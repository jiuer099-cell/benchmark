import importlib.util
import json
from pathlib import Path
import statistics
import sys
import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('bayestyper_stability', ROOT/'plugins/bayestyper/stability.py')
stability = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stability)


def noise_file(path, *, joint=False, separated=False):
    rows = ['Chain\tIteration\tHG002']
    retained=[]
    for c in range(1,21):
        for i in range(351):
            value = .2 + (c * .01 if separated else 0) + (.001 if i % 2 else -.001)
            rows.append(f'{c}\t{i}\t{value}')
            if i > 100: retained.append(value)
    if not joint: rows.append(f'0\t0\t{statistics.mean(retained)}')
    path.write_text('\n'.join(rows)+'\n', encoding='utf-8')


def test_stable_and_separated_chains_are_distinguished(tmp_path):
    p=tmp_path/'noise.txt'; noise_file(p)
    d=stability.analyze_file(p)
    assert d['status']=='NO_SCREEN_FLAG'
    assert d['samples']['HG002']['final_mean_absolute_difference']<1e-12
    noise_file(p,separated=True)
    assert stability.analyze_file(p)['status']=='REVIEW_REQUIRED'


def test_joint_format_does_not_need_sequential_final_row(tmp_path):
    p=tmp_path/'noise.txt'; noise_file(p,joint=True)
    assert stability.analyze_file(p)['status']=='PENDING'
    d=stability.analyze_file(p,noise_genotyping=True)
    assert d['status']=='NO_SCREEN_FLAG'
    assert d['samples']['HG002']['native_final_mean'] is None


def test_truncated_duplicate_and_nonfinite_chains_rejected(tmp_path):
    p=tmp_path/'noise.txt'; noise_file(p)
    original=p.read_text()
    p.write_text(original.replace('20\t350\t0.199\n',''))
    assert stability.analyze_file(p)['status']=='PENDING'
    p.write_text(original+'1\t1\t0.2\n')
    with pytest.raises(ValueError,match='Duplicate'):stability.analyze_file(p)
    p.write_text(original.replace('1\t1\t0.201','1\t1\tnan'))
    with pytest.raises(ValueError,match='Nonfinite'):stability.analyze_file(p)


def historical_run(tmp_path):
    if sys.platform == 'win32':
        tmp_path = Path('\\\\?\\' + str(tmp_path))
    base=tmp_path/'acceptance';base.mkdir()
    manifests=base/'manifests';manifests.mkdir()
    config=manifests/'config.yaml';config.write_text('frozen')
    root=base/'results/R/HG002/end_to_end_from_reads/bayestyper'
    (root/'raw').mkdir(parents=True);(root/'meta').mkdir()
    work=root/'.pgbench_attempts/abc/output/bayestyper-work'
    unit=work/'bayestyper_unit_1';unit.mkdir(parents=True)
    (root/'.pgbench_attempts/def').mkdir()
    native=work/'native.vcf.gz';native.write_bytes(b'native')
    (unit/'variant_clusters.bin').write_bytes(b'unit')
    noise_file(unit/'calls_noise_parameters.txt',separated=True)
    raw=root/'raw/calls.vcf';raw.write_bytes(b'raw')
    source=tmp_path/'code/plugins/bayestyper/run.py';source.parent.mkdir(parents=True);source.write_text('executed')
    frozen_hash=stability.sha(source)
    archive=base/'stability-fix/executed-source-before-stability-fix';archive.mkdir(parents=True)
    (archive/'run.py').write_text('executed')
    source.write_text('latest')
    (base/'run-spec.json').write_text(json.dumps({'run_id':'R','config':str(config),'configuration_sha256':stability.sha(config),'result':str(base/'results/R')}))
    (manifests/'source-state.json').write_text(json.dumps({'implementation_hashes':{source.as_posix():frozen_hash,'/unrelated/tool.py':'not-current'}}))
    attempt={'attempt_id':'abc','run_id':'R','status':'success','exit_code':0,'sample_id':'HG002','runner_path':source.as_posix(),'runner_sha256':frozen_hash,'output':{'sha256':stability.sha(raw)}}
    (root/'meta/attempt.abc.json').write_text(json.dumps(attempt))
    (root/'meta/attempt.def.json').write_text(json.dumps(dict(attempt,attempt_id='def',status='failed',exit_code=1)))
    audit={'output_sha256':stability.sha(raw),'artifacts':{'native_vcf':{'path':str(native),'sha256':stability.sha(native)}}}
    raw.with_suffix('.native-audit.json').write_text(json.dumps(audit))
    return base,root,source,archive


def test_historical_attempt_and_preserved_source_bindings(tmp_path):
    base,root,source,archive=historical_run(tmp_path)
    d=stability.inspect_run(base,'abc')
    assert d['status']=='REVIEW_REQUIRED' and d['attempt_id']=='abc'
    assert d['noise_bound_at_execution'] is False
    with pytest.raises(ValueError,match='successful'):stability.inspect_run(base,'def')
    (archive/'run.py').write_text('changed')
    with pytest.raises(ValueError,match='not preserved'):stability.inspect_run(base,'abc')


def test_historical_output_and_noise_tampering_detected(tmp_path):
    base,root,source,archive=historical_run(tmp_path)
    raw=root/'raw/calls.vcf';auditp=raw.with_suffix('.native-audit.json')
    audit=json.loads(auditp.read_text())
    noise=root/'.pgbench_attempts/abc/output/bayestyper-work/bayestyper_unit_1/calls_noise_parameters.txt'
    audit['artifacts']['noise_parameters_bayestyper_unit_1']={'path':str(noise),'sha256':stability.sha(noise)}
    auditp.write_text(json.dumps(audit))
    assert stability.inspect_run(base,'abc')['noise_bound_at_execution'] is True
    noise_file(noise)
    with pytest.raises(ValueError,match='Noise hash'):stability.inspect_run(base,'abc')
    raw.write_text('changed')
    with pytest.raises(ValueError,match='selected attempt'):stability.inspect_run(base,'abc')


def test_no_call_statistics_and_joint_command(tmp_path):
    sys.path.insert(0,str(ROOT/'plugins/bayestyper'))
    spec=importlib.util.spec_from_file_location('bayes_runner_screen_test',ROOT/'plugins/bayestyper/run.py')
    runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    p=tmp_path/'calls.vcf'
    p.write_text('#CHROM\nchr1\t1\tx\tA\tAA\t.\tAN0\t.\tGT\t./.\nchr1\t2\ty\tA\tAA\t.\tPASS\t.\tGT\t0/0\n')
    assert runner.genotype_summary(p)=={'total':2,'called':1,'no_call':1}
    default=runner.genotype_command('bayesTyper',Path('unit'),tmp_path,[],Path('calls'))
    joint=runner.genotype_command('bayesTyper',Path('unit'),tmp_path,[],Path('calls'),noise_genotyping=True)
    assert '--noise-genotyping' not in default
    assert joint==default+['--noise-genotyping']


@pytest.mark.parametrize('joint',[False,True])
def test_adapter_records_noise_artifacts_and_native_no_calls(tmp_path,monkeypatch,joint):
    import gzip
    sys.path.insert(0,str(ROOT/'plugins/bayestyper'))
    spec=importlib.util.spec_from_file_location('bayes_runner_audit_test',ROOT/'plugins/bayestyper/run.py')
    runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    reference=tmp_path/'ref.fa';reference.write_text('>chr1\nAAAAAAAAAA\n')
    fai=Path(str(reference)+'.fai');fai.write_text('chr1\t10\t6\t10\t11\n')
    panel=tmp_path/'panel.vcf';panel.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t3\tC1\tA\tAA\t.\t.\t.\n')
    assets=tmp_path/'assets';assets.mkdir()
    for name in ['canon.fa','decoy.fa','bayesTyper','bayesTyperTools','kmc']:(assets/name).write_text('asset')
    output_dir=tmp_path/'out';output_dir.mkdir()
    output=output_dir/'calls.vcf'
    inputs={'short_fastq_r1':tmp_path/'r1','short_fastq_r2':tmp_path/'r2','candidate_panel':panel,'reference':reference,'reference_index':fai,'tool_index':assets,'adapter_asset.binaries':assets}
    monkeypatch.setattr(runner,'declared_inputs',lambda:inputs)
    for k,val in {'PGBENCH_SAMPLE_ID':'HG002','PGBENCH_THREADS':'2','PGBENCH_MEMORY_MB':'3072','PGBENCH_RANDOM_SEED':'20260716','PGBENCH_OUTPUT_DIR':str(output_dir),'PGBENCH_OUTPUT_VCF':str(output)}.items():monkeypatch.setenv(k,val)
    commands=[]
    def tool(command,*,cwd=None):
        commands.append(command)
        if 'cluster' in command:
            unit=cwd/'bayestyper_unit_1';unit.mkdir();(unit/'variant_clusters.bin').write_bytes(b'unit')
        if 'genotype' in command:
            assert ('--noise-genotyping' in command)==joint
            prefix=Path(command[command.index('-o')+1])
            noise_file(prefix.parent/'calls_noise_parameters.txt',joint=joint,separated=True)
        if 'concat' in command:
            path=Path(command[command.index('-o')+1])
            with gzip.open(path,'wt') as f:f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\nchr1\t3\tC1\tA\tAA\t.\tAN0\t.\tGT\t./.\n')
    monkeypatch.setattr(runner,'run',tool)
    assert runner.main(noise_genotyping=joint)==0
    audit=json.loads(output.with_suffix('.native-audit.json').read_text())
    assert audit['genotype_summary']=={'total':1,'called':0,'no_call':1}
    assert audit['noise_stability']=='REVIEW_REQUIRED'
    assert audit['parameter_policy']['noise_genotyping']==joint
    for record in audit['artifacts'].values():assert stability.sha(record['path'])==record['sha256']
    assert 'noise_parameters_bayestyper_unit_1' in audit['artifacts']
    assert output.read_text().rstrip().endswith('GT\t./.')
