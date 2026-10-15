#!/usr/bin/env python3
"""Immutable shared BAM -> official RG conversion -> official Manta/GSE/GATK-SV."""
from __future__ import annotations
import argparse,datetime,gzip,json,os,pathlib,shutil,subprocess,sys
import run as base
import rg_compat

ROOT=pathlib.Path(__file__).resolve().parent


def save(path,value):
    path=pathlib.Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def status(out,phase,**extra):
    frozen=ROOT/'source.json'
    if frozen.is_file():
        for path,expected in json.loads(frozen.read_text())['sha256'].items():
            if base.digest(path)!=expected:raise ValueError('Frozen GATK implementation/config changed: '+path)
    row=dict(tool='gatk-sv',phase=phase,pid=os.getpid(),time_beijing=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),**extra)
    save(out/'run-status.json',row);print(json.dumps(row),flush=True)


def validate_inputs(config):
    base.check_contract(config)
    if config.get('compatibility_policy')!='official_RG_fill_and_original_reference_Manta_v1':raise ValueError('Missing explicit compatibility policy')
    if config.get('require_frozen_source_manifest') and not (ROOT/'source.json').is_file():
        raise ValueError('Missing production source freeze')
    for key in ('bam_or_cram','alignment_index','alignment_reference_dict','alignment_reference_fasta','alignment_reference_index',
                'alignment_freeze_manifest','runtime_root','download_root','helper_python','output_dir','samtools'):
        p=pathlib.Path(config[key]);config[key]=str(p.resolve())
        if any(c in config[key] for c in '\n\r\"$`,: ') :raise ValueError('Unsafe bind path: '+key)
        if key!='output_dir' and not p.exists():raise ValueError('Missing input: '+key)
    if pathlib.Path(config['bam_or_cram']).suffix!='.bam':raise ValueError('This compatibility policy requires BAM, not CRAM')
    if not 1<=config['compat_threads']<=config['threads']:raise ValueError('Invalid compatibility thread budget')
    out=pathlib.Path(config['output_dir'])
    for ro in ('runtime_root','download_root'):
        if out.is_relative_to(config[ro]):raise ValueError('Output overlaps immutable resources')
    return out


def stage_config(config,stage,alignment,index):
    x={**config,'output_dir':str(stage),'bam_or_cram':str(alignment),'alignment_index':str(index)}
    return x


def sandbox(config,command,extra_ro):
    result=base.sandbox_command(config,command)
    offset=result.index('--bind')
    additions=[]
    for path in sorted(set(str(pathlib.Path(p).resolve()) for p in extra_ro)):
        additions+=['--ro-bind',path,path]
    result[offset:offset]=additions
    return result


def run_stage(config,stage,name,workflow,inputs,alignment,index):
    stage.mkdir()
    for d in ('cache','tmp','home','cromwell-executions'):(stage/d).mkdir()
    shutil.copyfile(ROOT/'container_task.py',stage/'container_task.py')
    x=stage_config(config,stage,alignment,index)
    extra=[pathlib.Path(config['alignment_reference_fasta']).parent,pathlib.Path(alignment).parent,stage.parent]
    extra += [pathlib.Path(p) for p in config.get('resume_readonly_paths', [])]
    ctx={**x,'assets_dir':config['download_root']+'/assets','alignment_files':[str(alignment),str(index)],
         'extra_readonly_paths':[str(p) for p in extra],'task_memory_mb':config['memory_mb']-2048}
    save(stage/'task-context.json',ctx);save(stage/'inputs.json',inputs)
    (stage/'cromwell.conf').write_text(base.backend_config(stage,config['helper_python']))
    runtime=pathlib.Path(config['runtime_root']);download=pathlib.Path(config['download_root'])
    validation=subprocess.run([str(runtime/'bin/womtool'),'validate','-i',str(stage/'inputs.json'),str(workflow)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (stage/'wdl-validation.log').write_text(validation.stdout)
    if validation.returncode:raise ValueError(name+' WDL/input validation failed')
    command=[str(download/'runtime-downloads/java17/bin/java'),f'-Dconfig.file={stage}/cromwell.conf','-Xmx2g','-jar',
             str(download/'runtime-downloads/cromwell-91.jar'),'run',str(workflow),'-i',str(stage/'inputs.json'),'-m',str(stage/'metadata.json')]
    save(stage/'stage-plan.json',dict(workflow=str(workflow),workflow_sha256=base.digest(workflow),inputs_sha256=base.digest(stage/'inputs.json'),
        adapter_sha256=base.digest(__file__),backend_adapter_sha256=base.digest(base.__file__),
        container_task_sha256=base.digest(stage/'container_task.py'),backend_sha256=base.digest(stage/'cromwell.conf'),
        input_localization='copy',
        native_filters='unchanged',command=command))
    env={'PATH':'/usr/bin:/bin','HOME':str(stage/'home'),'TMPDIR':str(stage/'tmp'),'LANG':'C.UTF-8',
         'APPTAINER_CACHEDIR':str(stage/'cache'),'APPTAINER_TMPDIR':str(stage/'tmp'),'PYTHONDONTWRITEBYTECODE':'1'}
    base.execute(sandbox(x,command,extra),stage,config['timeout_seconds'],env)
    meta=json.loads((stage/'metadata.json').read_text())
    if meta.get('status')!='Succeeded':raise ValueError(name+' did not succeed')
    return meta


def input_schema(config,workflow):
    result=subprocess.check_output([config['runtime_root']+'/bin/womtool','inputs',str(workflow)],text=True)
    return json.loads(result)


def map_gse_inputs(primary,schema,manta,alignment,index):
    prefix='GATKSVPipelineSingleSample.'
    values={'GatherSampleEvidence.'+k[len(prefix):]:v for k,v in primary.items() if k.startswith(prefix) and 'GatherSampleEvidence.'+k[len(prefix):] in schema}
    values.update({'GatherSampleEvidence.bam_or_cram_file':str(alignment),'GatherSampleEvidence.bam_or_cram_index':str(index),
        'GatherSampleEvidence.manta_docker':None,'GatherSampleEvidence.melt_docker':None,
        'GatherSampleEvidence.manta_vcf_input':manta['Manta.vcf'],'GatherSampleEvidence.manta_vcf_index_input':manta['Manta.index'],
        'GatherSampleEvidence.run_module_metrics':primary.get(prefix+'run_sampleevidence_metrics',False),
        'GatherSampleEvidence.run_localize_reads':False})
    # Official shared-filesystem I/O switch; no read/variant filter changes.
    return values


def iter_calls(meta):
    for name,rows in meta.get('calls',{}).items():
        for row in rows:
            yield name,row
            if 'subWorkflowMetadata' in row:yield from iter_calls(row['subWorkflowMetadata'])


def evidence_gate(meta):
    outputs=meta['outputs'];proof={}
    for key in ('pesr_disc','pesr_split','pesr_sd'):
        p=pathlib.Path(outputs['GatherSampleEvidence.'+key])
        with gzip.open(p,'rt') as stream:
            if not stream.readline():raise ValueError('Empty '+key+' evidence; refuse downstream genotyping')
        proof[key]=dict(path=str(p),bytes=p.stat().st_size,sha256=base.digest(p))
    counts=pathlib.Path(outputs['GatherSampleEvidence.coverage_counts']);positive=False
    counts_opener = gzip.open if counts.suffix == '.gz' else open
    with counts_opener(counts, 'rt', encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('@') or line.startswith('CONTIG') or not line.strip():continue
            if int(line.rstrip().split('\t')[-1])>0:positive=True;break
    if not positive:raise ValueError('All-zero coverage counts; refuse downstream genotyping')
    proof['counts']=dict(path=str(counts),bytes=counts.stat().st_size,sha256=base.digest(counts),positive_counts=True)
    for name,row in iter_calls(meta):
        if name.endswith('RunCollectSVEvidence') or name.endswith('.CollectCounts'):
            text=pathlib.Path(row['stderr']).read_text()
            if 'Processed 0 total reads' in text:raise ValueError('Zero accepted reads in '+name)
    return dict(status='PASS',evidence=proof,default_read_filters_preserved=True)


def main(argv=None):
    if os.environ.get('PGBENCH_RESOLVED_INPUTS'):raise ValueError('Discovery is unsupported in fixed-panel PG-F1')
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=pathlib.Path,required=True)
    parser.add_argument('--execute',action='store_true');a=parser.parse_args(argv)
    config=json.loads(a.config.read_text());out=validate_inputs(config)
    original=pathlib.Path(config['bam_or_cram']);index=pathlib.Path(config['alignment_index'])
    header=subprocess.check_output([config['samtools'],'view','-H',str(original)],text=True)
    rg_compat.single_rg(header,config['sample_id'])
    reference_check=rg_compat.original_reference_check(header,pathlib.Path(config['alignment_reference_dict']).read_text(),base.parse_dictionary)
    primary,workflow=base.render_inputs(pathlib.Path(config['download_root']),config['sample_id'],original,index)
    broad_dict=pathlib.Path(primary['GATKSVPipelineSingleSample.reference_dict']).read_text()
    base.compatible_header(header,broad_dict,config['sample_id'],alignment_reference_dict=pathlib.Path(config['alignment_reference_dict']).read_text())
    subprocess.run([config['samtools'],'quickcheck','-v',str(original)],check=True)
    subprocess.run([config['samtools'],'view','-c','-X',str(original),str(index),'chr1:1-1'],stdout=subprocess.DEVNULL,check=True)
    if rg_compat.sha(config['alignment_freeze_manifest'])!=config['alignment_freeze_manifest_sha256']:
        raise ValueError('Source freeze manifest changed')
    freeze=json.loads(pathlib.Path(config['alignment_freeze_manifest']).read_text())
    if freeze['invariants']['alignment_records_sha256']!=config['source_alignment_records_sha256']:
        raise ValueError('Source record baseline changed')
    if not a.execute:
        print(json.dumps(dict(status='PREFLIGHT_PASS',source_unmodified=True,manta_reference=reference_check)));return 0
    out.mkdir(parents=True,exist_ok=False)
    plan=dict(config=config,source_identity=rg_compat.identity(original),source_index_identity=rg_compat.identity(index),
        adapter_sha256=base.digest(__file__),native_filters='unchanged',truth_used=False,primary_fixed_panel_score_eligible=False,
        compatibility_policy=config['compatibility_policy'],reference_check=reference_check,status='running')
    save(out/'attempt.json',plan)
    try:
        status(out,'rg_compatibility_conversion_and_full_audit')
        with (out/'compatibility.console.log').open('w') as log:
            alignment,new_index,audit=rg_compat.convert(config,out/'compatibility',header,log,ROOT/'rg_record_audit')
        plan['rg_compatibility_audit']=audit;save(out/'attempt.json',plan)
        primary,workflow=base.render_inputs(pathlib.Path(config['download_root']),config['sample_id'],alignment,new_index)
        source=workflow.parent
        prefix='GATKSVPipelineSingleSample.'
        manta={'Manta.bam_or_cram_file':str(alignment),'Manta.bam_or_cram_index':str(new_index),
            'Manta.sample_id':config['sample_id'],'Manta.reference_fasta':config['alignment_reference_fasta'],
            'Manta.reference_index':config['alignment_reference_index'],'Manta.region_bed':primary[prefix+'manta_region_bed'],
            'Manta.region_bed_index':primary[prefix+'manta_region_bed_index'],'Manta.manta_docker':primary[prefix+'manta_docker']}
        for key in ('jobs_per_cpu','mem_gb_per_job'):
            if prefix+'manta_'+key in primary:manta['Manta.'+key]=primary[prefix+'manta_'+key]
        if prefix+'runtime_attr_manta' in primary:manta['Manta.runtime_attr_override']=primary[prefix+'runtime_attr_manta']
        status(out,'manta_original_alignment_reference')
        mm=run_stage(config,out/'manta','Manta',source/'Manta.wdl',manta,alignment,new_index)
        schema=input_schema(config,source/'GatherSampleEvidence.wdl')
        gse=map_gse_inputs(primary,schema,mm['outputs'],alignment,new_index)
        status(out,'gather_sample_evidence_default_filters')
        gm=run_stage(config,out/'evidence','GatherSampleEvidence',source/'GatherSampleEvidence.wdl',gse,alignment,new_index)
        proof=evidence_gate(gm);save(out/'evidence-validation.json',proof)
        mapping={'case_manta_vcf':'manta_vcf','case_wham_vcf':'wham_vcf','case_scramble_vcf':'scramble_vcf',
                 'case_counts_file':'coverage_counts','case_pe_file':'pesr_disc','case_sr_file':'pesr_split','case_sd_file':'pesr_sd'}
        for target,key in mapping.items():
            value=gm['outputs'].get('GatherSampleEvidence.'+key)
            if not value:raise ValueError('Missing official precomputed evidence: '+key)
            primary[prefix+target]=value
        status(out,'official_single_sample_downstream')
        pm=run_stage(config,out/'pipeline','GATKSVPipelineSingleSample',workflow,primary,alignment,new_index)
        plan['native_output']=base.publish_native(pm,out/'pipeline',config['sample_id']);plan['status']='succeeded'
        status(out,'complete',native_output=plan['native_output'])
    except BaseException as exc:
        plan['status']='failed';plan['error']=str(exc);status(out,'failed',error=str(exc));raise
    finally:
        plan['finished_at_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();save(out/'attempt.json',plan)
    return 0


if __name__=='__main__':raise SystemExit(main())
