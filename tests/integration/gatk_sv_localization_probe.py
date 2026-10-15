"""Run real Cromwell/Apptainer input-localization regression on prepared runtime.

Uses a completed counts file read-only; writes only a new probe directory.
"""
import argparse
import gzip
import json
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / 'plugins/gatk-sv'))
import run_rg_compatible as adapter

WDL = '''version 1.0
workflow LocalizationProbe {
  input {
    File counts
    String image
  }
  call Produce { input: counts=counts, image=image }
  call Consume { input: counts=Produce.produced, image=image }
  output { File uncompressed=Consume.uncompressed }
}
task Produce {
  input {
    File counts
    String image
  }
  command <<<
    set -euo pipefail
    test ! -L ~{counts}
    cp ~{counts} producer.tsv.gz
  >>>
  output { File produced="producer.tsv.gz" }
  runtime {
    docker: image
    cpu: 1
    memory_gb: 1.0
  }
}
task Consume {
  input {
    File counts
    String image
  }
  command <<<
    set -euo pipefail
    test ! -L ~{counts}
    printf '%s\\n' '~{counts}' | xargs -l1 -P0 gunzip
    cp ~{sub(counts, "\\\\.gz$", "")} uncompressed.tsv
  >>>
  output { File uncompressed="uncompressed.tsv" }
  runtime {
    docker: image
    cpu: 1
    memory_gb: 1.0
  }
}
'''


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--config',required=True,type=pathlib.Path)
    p.add_argument('--counts',required=True,type=pathlib.Path)
    p.add_argument('--output',required=True,type=pathlib.Path)
    a=p.parse_args()
    config=json.loads(a.config.read_text())
    original=a.counts.resolve(strict=True)
    before=adapter.base.digest(original)
    a.output.mkdir(parents=True,exist_ok=False)
    fixture=a.output/'counts.tsv.gz'
    fixture.write_bytes(original.read_bytes())
    config['output_dir']=str(a.output)
    config['resume_readonly_paths']=[]
    image='us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT'
    records=json.loads((pathlib.Path(config['runtime_root'])/'sif_images.json').read_text())['images']
    sif=next(r['sif'] for r in records if image in (r['image'],r['actual_image']))
    link=a.output/'linked.tsv.gz';link.symlink_to(fixture)
    cmd=[config['runtime_root']+'/apptainer/bin/apptainer','exec','--userns','--cleanenv','--containall','--no-home',
         '--pwd',str(a.output),'--bind',str(a.output)+':'+str(a.output),sif,'gunzip',str(link)]
    bad=subprocess.run(cmd,capture_output=True,text=True)
    assert bad.returncode!=0 and 'Too many levels of symbolic links' in bad.stderr,bad.stderr
    assert adapter.base.digest(fixture)==before
    wdl=a.output/'probe.wdl';wdl.write_text(WDL)
    meta=adapter.run_stage(config,a.output/'pipeline','LocalizationProbe',wdl,
        {'LocalizationProbe.counts':str(fixture),'LocalizationProbe.image':image},
        config['bam_or_cram'],config['alignment_index'])
    with gzip.open(original,'rb') as f: expected=f.read()
    assert pathlib.Path(meta['outputs']['LocalizationProbe.uncompressed']).read_bytes()==expected
    assert adapter.base.digest(original)==before
    assert adapter.base.digest(fixture)==before
    producer=list((a.output/'pipeline/cromwell-executions').rglob('producer.tsv.gz'))
    assert any('/execution/' in str(f) and f.is_file() and not f.is_symlink() for f in producer)
    report=dict(status='PASS',workflow_status=meta['status'],old_soft_link_returncode=bad.returncode,
                old_soft_link_stderr=bad.stderr,source_sha256=before,source_unchanged=True,
                consumer_output_matches_decompressed_source=True,backend_sha256=adapter.base.digest(a.output/'pipeline/cromwell.conf'))
    adapter.save(a.output/'regression-report.json',report)
    print(json.dumps(report))


if __name__=='__main__': main()
