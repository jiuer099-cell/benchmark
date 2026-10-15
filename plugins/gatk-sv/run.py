#!/usr/bin/env python3
"""Standalone GATK-SV discovery experiment; never a fixed-panel PG-F1 plugin."""
from __future__ import annotations
import argparse,datetime,gzip,hashlib,json,os,pathlib,re,shutil,signal,subprocess,sys

ROOT=pathlib.Path(__file__).resolve().parent
INTERFACE='gatk_sv_single_sample_v1'


def digest(path):
    h=hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def check_contract(config):
    if config.get('experiment_interface')!=INTERFACE:raise ValueError('Unknown experiment interface')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',config['sample_id']):raise ValueError('Invalid sample ID')
    for field in ('threads','memory_mb','timeout_seconds'):
        if type(config[field]) is not int or config[field]<1:raise ValueError('Invalid resource limit: '+field)
    if config['memory_mb']<4096:raise ValueError('At least 4096 MB is required for the engine and tasks')
    forbidden={'truth','truth_vcf','candidate_panel','panel_vcf','primary_score','PG-F1'}
    if forbidden.intersection(config):raise ValueError('Discovery config cannot contain truth, canonical-panel projection or primary scoring inputs')


def parse_dictionary(text):
    result={}
    for line in text.splitlines():
        if line.startswith('@SQ\t'):
            fields=dict(f.split(':',1) for f in line.split('\t')[1:] if ':' in f)
            if fields['SN'] in result:raise ValueError('Duplicate alignment/reference contig')
            result[fields['SN']]=fields
    if not result:raise ValueError('Missing @SQ sequence dictionary')
    return result


def compatible_header(header, reference_dict, sample, *, cram=False, alignment_reference_dict=None):
    expected=parse_dictionary(reference_dict);observed=parse_dictionary(header)
    alignment_reference=parse_dictionary(alignment_reference_dict) if alignment_reference_dict else {}
    for name,row in observed.items():
        ref=expected.get(name)
        if ref is None or row['LN']!=ref['LN']:raise ValueError('Alignment/reference contig mismatch: '+name)
        if cram and not row.get('M5'):raise ValueError('CRAM requires @SQ M5 reference checksums')
        if not row.get('M5'):
            declared=alignment_reference.get(name,{})
            if declared.get('LN')!=row['LN'] or not declared.get('M5'):
                raise ValueError('BAM lacks @SQ M5; provide the frozen alignment_reference_dict: '+name)
            row['M5']=declared['M5']
        if row.get('M5') and row['M5'].lower()!=ref.get('M5','').lower():raise ValueError('Alignment/reference M5 mismatch: '+name)
    samples={f[3:] for line in header.splitlines() if line.startswith('@RG\t') for f in line.split('\t')[1:] if f.startswith('SM:')}
    if samples!={sample}:raise ValueError('Alignment read-group sample IDs do not match the requested single sample')
    return {'contigs_checked':len(observed),'sample_id':sample,'cram_md5_required':cram}


def check_alignment_files(container, alignment, index, reference, header):
    """Reject truncated inputs and unusable explicit indexes before execution.

    Only a one-base interval is read; idxstats could silently fall back to
    scanning the whole alignment when the index is absent or invalid.
    """
    contig=next(iter(parse_dictionary(header)))
    commands=[[*container,'samtools','quickcheck','-v',str(alignment)],
              [*container,'samtools','view','-c','-X','-T',str(reference),
               str(alignment),str(index),contig+':1-1']]
    for command in commands:
        result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=60)
        if result.returncode:
            raise ValueError('Alignment/index preflight failed: '+result.stderr.strip())
    return {'quickcheck':True,'explicit_index_query':True,'query_interval':contig+':1-1',
            'scope':'Header/EOF and bounded indexed read; not a complete alignment/index content audit'}


def render_inputs(download, sample, alignment, index):
    from jinja2 import Environment,StrictUndefined
    source=download/'source/broadinstitute-gatk-sv-8b41fa8'
    context=json.loads((download/'assets/selected_inputs.local.json').read_text())
    context['single_sample']=dict(sample_id=sample,name=sample,bam_or_cram_file=str(alignment))
    template=(source/'inputs/templates/test/GATKSVPipelineSingleSample/GATKSVPipelineSingleSample.json.tmpl').read_text()
    inputs=json.loads(Environment(undefined=StrictUndefined).from_string(template).render(**context))
    prefix='GATKSVPipelineSingleSample.'
    inputs[prefix+'bam_or_cram_index']=str(index)
    inputs[prefix+'use_melt']=False;inputs[prefix+'run_vcf_qc']=False
    return inputs,source/'wdl/GATKSVPipelineSingleSample.wdl'


def backend_config(output, python):
    # No global runtime/cromwell execution directory or cache is reused.
    # Official gCNV tasks gunzip localized inputs in place. Soft links are
    # refused by gunzip (ELOOP); hard links can also expose upstream data to
    # writes. Give every call its own regular file, including retries.
    task=output/'container_task.py';ctx=output/'task-context.json'
    return f'''include required(classpath("application"))
docker.hash-lookup.enabled = false
call-caching.enabled = false
backend {{
 default = "Apptainer"
 providers.Apptainer {{
  actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
  config {{
   root = "{output}/cromwell-executions"
   concurrent-job-limit = 1
   run-in-background = true
   runtime-attributes = """
     String? docker
     Int cpu = 2
     Float memory_gb = 4.0
     String? disks
     Int? bootDiskSizeGb
     Int? preemptible
     Int? maxRetries
   """
   submit = "/bin/bash ${{script}}"
   submit-docker = """
     {python} {task} --context {ctx} --image "${{docker}}" --cwd "${{cwd}}" --container-cwd "${{docker_cwd}}" --shell "${{job_shell}}" --script "${{docker_script}}" --pid-file "${{docker_cid}}" --cpu "${{cpu}}" --memory-gb "${{memory_gb}}"
   """
   kill-docker = "/bin/kill -TERM -- -$(cat ${{docker_cid}})"
   filesystems.local {{
    localization = ["copy"]
    caching.duplication-strategy = ["copy"]
    docker.allow-soft-links = false
   }}
  }}
 }}
}}
'''


def sandbox_command(config, command):
    runtime=pathlib.Path(config['runtime_root']);download=pathlib.Path(config['download_root'])
    python=pathlib.Path(config['helper_python']);out=pathlib.Path(config['output_dir'])
    ro=[pathlib.Path(p) for p in ['/usr','/bin','/lib','/lib64','/sys','/etc/passwd','/etc/group','/etc/hosts','/etc/nsswitch.conf','/etc/localtime'] if pathlib.Path(p).exists()]
    ro += [runtime/'apptainer', runtime/'images', runtime/'sif_images.json', download/'runtime-downloads',
           download/'source/broadinstitute-gatk-sv-8b41fa8',download/'assets',python.resolve().parent.parent,
           pathlib.Path(config['bam_or_cram']),pathlib.Path(config['alignment_index'])]
    # bwrap creates bind parent directories itself. No broad home/data mount.
    cmd=['bwrap','--die-with-parent','--unshare-all','--unshare-net','--new-session','--proc','/proc','--dev','/dev','--dev-bind','/dev/fuse','/dev/fuse']
    for path in ro:cmd+=['--ro-bind',str(path),str(path)]
    cmd+=['--bind',str(out),str(out),'--chdir',str(out),'--',*command]
    return cmd


def execute(command,output,timeout,env):
    with (output/'cromwell.log').open('w') as log:
        p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        try:
            rc=p.wait(timeout=timeout)
            if rc:raise RuntimeError(f'Cromwell exited with status {rc}; see cromwell.log')
        except BaseException:
            if p.poll() is None:
                os.killpg(p.pid,signal.SIGTERM)
                try:p.wait(timeout=15)
                except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
            raise


def publish_native(metadata,output,sample):
    if metadata.get('status')!='Succeeded':raise ValueError('Workflow did not succeed; no native output is published')
    source=pathlib.Path(metadata['outputs']['GATKSVPipelineSingleSample.final_vcf']).resolve()
    source.relative_to(output.resolve())
    opener=gzip.open if source.suffix=='.gz' else open
    with opener(source,'rt') as stream:
        for line in stream:
            if line.startswith('#CHROM'):
                if sample not in line.rstrip().split('\t')[9:]:raise ValueError('Final native VCF lacks the requested sample')
                break
        else:raise ValueError('Final native VCF has no sample header')
    destination=output/'raw'/('calls.vcf.gz' if source.suffix=='.gz' else 'calls.vcf')
    destination.parent.mkdir(exist_ok=True)
    if destination.exists():raise ValueError('Refusing to overwrite published native output')
    shutil.copyfile(source,destination)
    index=metadata['outputs'].get('GATKSVPipelineSingleSample.final_vcf_idx')
    if index:
        path=pathlib.Path(index).resolve();path.relative_to(output.resolve())
        if path.suffix not in ('.csi','.tbi','.idx'):raise ValueError('Unexpected native VCF index extension')
        shutil.copyfile(path,pathlib.Path(str(destination)+path.suffix))
    return {'path':str(destination),'sha256':digest(destination),'source_path':str(source),
            'native_records_and_filters_preserved':True,'canonical_projection_applied':False}


def main(argv=None):
    if os.environ.get('PGBENCH_RESOLVED_INPUTS'):
        raise ValueError('GATK-SV discovery adapter is unsupported in the fixed-panel PG-F1 executor')
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=pathlib.Path,required=True)
    mode=p.add_mutually_exclusive_group();mode.add_argument('--prepare-only',action='store_true');mode.add_argument('--execute',action='store_true')
    a=p.parse_args(argv);config=json.loads(a.config.read_text());check_contract(config)
    for key in ('runtime_root','download_root','helper_python','bam_or_cram','alignment_index','output_dir'):
        config[key]=str(pathlib.Path(config[key]).resolve())
        if any(c in config[key] for c in '\n\r\"$`,:') or ' ' in config[key]:raise ValueError('Runtime/experiment paths must contain no shell/bind metacharacters or spaces')
    runtime=pathlib.Path(config['runtime_root']);download=pathlib.Path(config['download_root']);out=pathlib.Path(config['output_dir'])
    alignment=pathlib.Path(config['bam_or_cram']);index=pathlib.Path(config['alignment_index'])
    for readonly_root in (runtime,download):
        if out.is_relative_to(readonly_root):raise ValueError('Output must be outside immutable runtime/resources')
    for file in (alignment,index):
        if not file.is_file() or not file.stat().st_size:raise ValueError('Missing/empty alignment input: '+str(file))
    if alignment.suffix not in ('.bam','.cram'):raise ValueError('Expected a BAM or CRAM file')
    if not shutil.which('bwrap'):raise ValueError('bwrap is required for network-isolated discovery execution')
    verified=json.loads((runtime/'runtime_verification.json').read_text())
    if not verified['success'] or verified['versions']['cromwell']!='cromwell 91':raise ValueError('The tested Cromwell91/Apptainer runtime is required')
    inputs,workflow=render_inputs(download,config['sample_id'],alignment,index)
    resource=inputs['GATKSVPipelineSingleSample.reference_dict']
    records=json.loads((runtime/'sif_images.json').read_text())['images']
    mini=next(r['sif'] for r in records if r['key']=='sv_base_mini_docker')
    container=[str(runtime/'bin/apptainer'),'exec','--userns','--cleanenv','--containall','--no-home','--pwd','/',
               '--bind',f'{alignment}:{alignment}:ro','--bind',f'{index}:{index}:ro',
               '--bind',f'{download}/assets:{download}/assets:ro',mini]
    reference=inputs['GATKSVPipelineSingleSample.reference_fasta']
    header_cmd=[*container,'samtools','view','-H','-T',reference,str(alignment)]
    header=subprocess.check_output(header_cmd,text=True,timeout=60)
    declared_reference=config.get('alignment_reference_dict')
    reference_check=compatible_header(header,pathlib.Path(resource).read_text(),config['sample_id'],cram=alignment.suffix=='.cram',
                                      alignment_reference_dict=pathlib.Path(declared_reference).read_text() if declared_reference else None)
    alignment_check=check_alignment_files(container,alignment,index,reference,header)
    out.mkdir(parents=True,exist_ok=False)
    for d in ('cache','tmp','home','cromwell-executions'):(out/d).mkdir()
    shutil.copyfile(ROOT/'container_task.py',out/'container_task.py')
    ctx={**config,'assets_dir':str(download/'assets'),'alignment_files':[str(alignment),str(index)],'task_memory_mb':config['memory_mb']-2048}
    (out/'task-context.json').write_text(json.dumps(ctx,indent=2)+'\n')
    (out/'inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    (out/'cromwell.conf').write_text(backend_config(out,config['helper_python']))
    validation=subprocess.run([str(runtime/'bin/womtool'),'validate','-i',str(out/'inputs.json'),str(workflow)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (out/'wdl-validation.log').write_text(validation.stdout)
    if validation.returncode:raise ValueError('WDL/input validation failed; see wdl-validation.log')
    java=download/'runtime-downloads/java17/bin/java';jar=download/'runtime-downloads/cromwell-91.jar'
    command=[str(java),f'-Dconfig.file={out}/cromwell.conf','-Xmx2g','-jar',str(jar),'run',str(workflow),'-i',str(out/'inputs.json'),'-m',str(out/'metadata.json')]
    plan={'experiment_interface':INTERFACE,'scope':'sv_discovery_and_genotyping','primary_fixed_panel_score_eligible':False,
          'truth_used':False,'alignment_header_compatibility':reference_check,'alignment_integrity':alignment_check,'config':config,'inputs_sha256':digest(out/'inputs.json'),
          'adapter_sha256':digest(__file__),'images_manifest_sha256':digest(runtime/'sif_images.json'),
          'container_task_sha256':digest(out/'container_task.py'),'backend_sha256':digest(out/'cromwell.conf'),
          'workflow_sha256':digest(workflow),
          'command':sandbox_command(config,command),'native_filters':'unchanged','status':'prepared',
          'memory_policy':'Reject task memory requests exceeding task budget; CPU affinity capped; no OS RSS cgroup enforcement.'}
    (out/'attempt.json').write_text(json.dumps(plan,indent=2)+'\n')
    print('Prepared isolated GATK-SV discovery attempt: '+str(out),flush=True)
    if not a.execute:return 0
    env={'PATH':'/usr/bin:/bin','HOME':str(out/'home'),'TMPDIR':str(out/'tmp'),'LANG':'C.UTF-8',
         'APPTAINER_CACHEDIR':str(out/'cache'),'APPTAINER_TMPDIR':str(out/'tmp'),'PYTHONDONTWRITEBYTECODE':'1'}
    plan['status']='running';(out/'attempt.json').write_text(json.dumps(plan,indent=2)+'\n')
    try:
        execute(plan['command'],out,config['timeout_seconds'],env)
        metadata=json.loads((out/'metadata.json').read_text())
        plan['native_output']=publish_native(metadata,out,config['sample_id']);plan['status']='succeeded'
    except BaseException as exc:plan['status']='failed';plan['error']=str(exc);raise
    finally:
        plan['finished_at_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        (out/'attempt.json').write_text(json.dumps(plan,indent=2)+'\n')
    return 0


if __name__=='__main__':raise SystemExit(main())
