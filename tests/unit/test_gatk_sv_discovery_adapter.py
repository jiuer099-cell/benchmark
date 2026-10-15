from __future__ import annotations
import importlib.util,json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[2]/'plugins/gatk-sv'
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,ROOT/file);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
runner=load('gatk_discovery','run.py');task=load('gatk_container','container_task.py')

def test_discovery_refuses_fixed_panel_executor(monkeypatch):
    monkeypatch.setenv('PGBENCH_RESOLVED_INPUTS','a.json')
    with pytest.raises(ValueError,match='unsupported'):runner.main([])

@pytest.mark.parametrize('observed,cram',[
    ('@SQ\tSN:chr1\tLN:11\tM5:a\n@RG\tID:r\tSM:HG002',False),
    ('@SQ\tSN:chr1\tLN:10\tM5:b\n@RG\tID:r\tSM:HG002',False),
    ('@SQ\tSN:chr1\tLN:10\n@RG\tID:r\tSM:HG002',True),
    ('@SQ\tSN:chr1\tLN:10\tM5:a\n@RG\tID:r\tSM:OTHER',False),
])
def test_alignment_reference_and_sample_mismatch_fails(observed,cram):
    with pytest.raises(ValueError):runner.compatible_header(observed,'@SQ\tSN:chr1\tLN:10\tM5:a','HG002',cram=cram)

def test_valid_alignment_header():
    header='@SQ\tSN:chr1\tLN:10\tM5:a\n@RG\tID:r\tSM:HG002'
    assert runner.compatible_header(header,header,'HG002',cram=True)['contigs_checked']==1

def test_bam_without_md5_requires_frozen_alignment_reference():
    header='@SQ\tSN:chr1\tLN:10\n@RG\tID:r\tSM:HG002'
    reference='@SQ\tSN:chr1\tLN:10\tM5:a'
    with pytest.raises(ValueError,match='alignment_reference_dict'):runner.compatible_header(header,reference,'HG002')
    assert runner.compatible_header(header,reference,'HG002',alignment_reference_dict=reference)['contigs_checked']==1

def test_image_mapping_only_allows_prepared_images():
    rows=[dict(image='old/ubuntu',actual_image='ubuntu@sha256:a',sif='/image.sif')]
    assert task.image_for(rows,'old/ubuntu')=='/image.sif'
    with pytest.raises(ValueError):task.image_for(rows,'unprepared/image')

def test_startup_failure_writes_rc_instead_of_hanging(tmp_path):
    cwd=tmp_path/'call';cwd.mkdir()
    ctx=tmp_path/'context.json';ctx.write_text(json.dumps(dict(output_dir=str(tmp_path),task_memory_mb=1000)))
    rc=task.main(['--context',str(ctx),'--image','unknown','--cwd',str(cwd),'--container-cwd','/call','--shell','/bin/bash','--script','/call/script','--pid-file',str(cwd/'pid'),'--cpu','2','--memory-gb','2'])
    assert rc==125 and (cwd/'execution/rc').read_text()=='125\n'

def test_failed_workflow_cannot_publish_native_output(tmp_path):
    with pytest.raises(ValueError,match='did not succeed'):runner.publish_native({'status':'Failed'},tmp_path,'HG002')

def test_native_publish_is_byte_preserving(tmp_path):
    native=tmp_path/'native.vcf';native.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\nchr1\t3\tSV1\tA\t<DEL>\t8\tLowQual\tSVTYPE=DEL\tGT\t0/1\n')
    metadata={'status':'Succeeded','outputs':{'GATKSVPipelineSingleSample.final_vcf':str(native)}}
    audit=runner.publish_native(metadata,tmp_path,'HG002')
    assert Path(audit['path']).read_bytes()==native.read_bytes()
    assert not audit['canonical_projection_applied']

def test_truth_config_rejected():
    x=dict(experiment_interface=runner.INTERFACE,sample_id='HG002',threads=2,memory_mb=4096,timeout_seconds=10,truth_vcf='secret')
    with pytest.raises(ValueError,match='truth'):runner.check_contract(x)


def test_index_preflight_uses_explicit_index_and_bounded_region(monkeypatch):
    commands=[]
    def run(command,**kwargs):
        commands.append(command)
        return runner.subprocess.CompletedProcess(command,0,stdout='0\n',stderr='')
    monkeypatch.setattr(runner.subprocess,'run',run)
    check=runner.check_alignment_files(['container'],'a.bam','custom.bai','ref.fa','@SQ\tSN:chr1\tLN:100')
    assert commands[0]==['container','samtools','quickcheck','-v','a.bam']
    assert commands[1]==['container','samtools','view','-c','-X','-T','ref.fa','a.bam','custom.bai','chr1:1-1']
    assert check['explicit_index_query']


@pytest.mark.parametrize('fail_on',[1,2])
def test_invalid_alignment_or_index_fails_preflight(monkeypatch,fail_on):
    calls=[]
    def run(command,**kwargs):
        calls.append(command)
        return runner.subprocess.CompletedProcess(command,int(len(calls)==fail_on),stdout='',stderr='invalid input')
    monkeypatch.setattr(runner.subprocess,'run',run)
    with pytest.raises(ValueError,match='preflight failed'):
        runner.check_alignment_files(['container'],'a.bam','custom.bai','ref.fa','@SQ\tSN:chr1\tLN:100')
    assert len(calls)==fail_on
