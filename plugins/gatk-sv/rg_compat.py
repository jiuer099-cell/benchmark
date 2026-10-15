"""Official samtools RG conversion; immutable source, full record-equivalence gate."""
from __future__ import annotations
import hashlib,json,pathlib,subprocess,time


def identity(path):
    s=pathlib.Path(path).stat()
    return dict(device=s.st_dev,inode=s.st_ino,size=s.st_size,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns)


def sha(path):
    h=hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def single_rg(header,sample):
    rows=[]
    for line in header.splitlines():
        if line.startswith('@RG\t'):
            rows.append(dict(t.split(':',1) for t in line.split('\t')[1:] if ':' in t))
    if len(rows)!=1 or rows[0].get('SM')!=sample or not rows[0].get('ID'):
        raise ValueError('RG fill requires exactly one declared source read group for the requested sample')
    return rows[0]['ID']


def original_reference_check(header,reference_dict,parse):
    observed=parse(header);expected=parse(reference_dict)
    if list(observed)!=list(expected):raise ValueError('Manta reference contig set/order must exactly match frozen BAM')
    for name,row in observed.items():
        if row['LN']!=expected[name]['LN']:raise ValueError('Manta reference length mismatch: '+name)
        if row.get('M5') and row['M5'].lower()!=expected[name].get('M5','').lower():
            raise ValueError('Manta reference sequence mismatch: '+name)
    return dict(contigs=len(observed),contig_order_identical=True)


def convert(config,directory,header,log,auditor):
    source=pathlib.Path(config['bam_or_cram']);directory=pathlib.Path(directory)
    directory.mkdir(parents=True,exist_ok=False)
    manifest=pathlib.Path(config['alignment_freeze_manifest'])
    if sha(manifest)!=config['alignment_freeze_manifest_sha256']:raise ValueError('Frozen source manifest hash changed')
    frozen=json.loads(manifest.read_text());record_hash=frozen['invariants']['alignment_records_sha256']
    if record_hash!=config['source_alignment_records_sha256']:raise ValueError('Frozen record baseline mismatch')
    before=identity(source);index_before=identity(config['alignment_index'])
    if before['size']!=frozen['outputs']['bam']['bytes']:raise ValueError('Frozen BAM size changed')
    rg=single_rg(header,config['sample_id']);samtools=config['samtools']
    version=subprocess.check_output([samtools,'--version'],text=True).splitlines()[0]
    if version!='samtools 1.10':raise ValueError('Record serialization baseline requires pinned samtools 1.10')
    dest=directory/(config['sample_id']+'.rg-compatible.bam');partial=dest.with_suffix('.bam.partial')
    command=[samtools,'addreplacerg','-@',str(config['compat_threads']),'-m','orphan_only','-R',rg,'-o',str(partial),'-O','BAM',str(source)]
    (directory/'conversion-plan.json').write_text(json.dumps(dict(source=str(source),source_identity=before,source_index_identity=index_before,
        source_sha256_declared=frozen['outputs']['bam']['sha256'],source_records_sha256=record_hash,command=command,
        samtools_sha256=sha(samtools),samtools_version=version,operation='fill_missing_RG_only',source_modified=False),indent=2)+'\n')
    started=time.time()
    subprocess.run(command,stdout=log,stderr=log,check=True)
    subprocess.run([samtools,'quickcheck','-v',str(partial)],stdout=log,stderr=log,check=True)
    read=subprocess.Popen([samtools,'view','-@',str(config['compat_threads']),str(partial)],stdout=subprocess.PIPE,stderr=log)
    try:
        audit=subprocess.run([str(auditor),rg],stdin=read.stdout,stdout=subprocess.PIPE,stderr=log,text=True,check=True)
    finally:
        read.stdout.close()
        read_rc=read.wait()
    if read_rc:raise ValueError('Complete converted BAM record traversal failed')
    result=json.loads(audit.stdout)
    if result['normalized_records_sha256']!=record_hash:raise ValueError('Non-RG records changed: refuse biological workflow')
    if identity(source)!=before or identity(config['alignment_index'])!=index_before:raise ValueError('Frozen source changed during compatibility conversion')
    # Same @SQ and original @RG; only the official samtools @PG line may be added.
    after_header=subprocess.check_output([samtools,'view','-H',str(partial)],text=True)
    fixed=lambda text:[l for l in text.splitlines() if not l.startswith('@PG\t')]
    if fixed(after_header)!=fixed(header):raise ValueError('Unexpected non-PG header changes')
    partial.rename(dest)
    subprocess.run([samtools,'index','-@',str(config['compat_threads']),str(dest)],stdout=log,stderr=log,check=True)
    result.update(status='PASS',source=str(source),source_identity_unchanged=True,source_index_identity_unchanged=True,
        source_sha256_declared=frozen['outputs']['bam']['sha256'],source_byte_sha256_rehashed_this_run=False,
        source_freeze_manifest_sha256=sha(manifest),source_records_sha256=record_hash,
        source_baseline_kind='pre-existing_frozen_full_SAM_records_sha256',
        output=str(dest),output_sha256=sha(dest),index_sha256=sha(str(dest)+'.bai'),
        header_change='official_samtools_PG_only',only_record_change='RG:Z:'+rg,
        every_alignment_field_and_non_RG_aux_preserved=True,elapsed_seconds=time.time()-started)
    (directory/'rg-compatibility-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    return dest,pathlib.Path(str(dest)+'.bai'),result
