#!/usr/bin/env python3
"""Run one pinned GATK-SV SIF with experiment-private writable paths."""
from __future__ import annotations
import argparse,json,os,pathlib,signal,subprocess,sys


def image_for(records, image):
    matches=[r for r in records if image in (r['image'],r['actual_image'])]
    if not matches:raise ValueError('Image is not in the pinned offline manifest: '+image)
    return matches[0]['sif']


def main(argv=None):
    p=argparse.ArgumentParser()
    for key in ('context','image','cwd','container-cwd','shell','script','pid-file'):p.add_argument('--'+key,required=True)
    p.add_argument('--cpu',type=int,required=True);p.add_argument('--memory-gb',type=float,required=True)
    a=p.parse_args(argv);cwd=pathlib.Path(a.cwd);child=None;safe_cwd=False
    try:
        x=json.loads(pathlib.Path(a.context).read_text())
        cwd.resolve().relative_to(pathlib.Path(x['output_dir']))
        safe_cwd=True
        if a.memory_gb*1000 > x['task_memory_mb']:
            raise ValueError(f'Task requests {a.memory_gb} GB, exceeding the declared task budget; raise memory_mb before starting a new attempt.')
        records=json.loads(pathlib.Path(x['runtime_root'],'sif_images.json').read_text())['images']
        sif=image_for(records,a.image)
        cmd=[str(pathlib.Path(x['runtime_root'],'apptainer/bin/apptainer')),'exec','--userns','--cleanenv','--containall','--no-home','--env','CLOUDSDK_CONFIG=/tmp/gatk-sv-cloud-sdk','--pwd',a.container_cwd]
        binds=[f'{path}:{path}:ro' for path in x.get('extra_readonly_paths',[])]
        binds += [str(cwd)+':'+a.container_cwd, x['output_dir']+':'+x['output_dir'], x['assets_dir']+':'+x['assets_dir']+':ro']
        binds += [f'{path}:{path}:ro' for path in x['alignment_files']]
        for bind in binds:cmd+=['--bind',bind]
        cmd += [sif,a.shell,a.script]
        cpus=sorted(os.sched_getaffinity(0))[:max(1,min(a.cpu,x['threads']))]
        os.sched_setaffinity(0,cpus)
        child=subprocess.Popen(cmd,cwd=cwd,start_new_session=True)
        pathlib.Path(a.pid_file).write_text(str(child.pid)+'\n')
        def stop(signum,frame):
            if child.poll() is None:os.killpg(child.pid,signum)
        signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
        status=child.wait();rc=status if status>=0 else 128-status
    except Exception as exc:
        print(str(exc),file=sys.stderr);rc=125
    rcfile=cwd/'execution/rc'
    if rc and safe_cwd and not rcfile.exists():
        rcfile.parent.mkdir(parents=True,exist_ok=True);rcfile.write_text(str(rc)+'\n')
    return rc


if __name__=='__main__':raise SystemExit(main())
