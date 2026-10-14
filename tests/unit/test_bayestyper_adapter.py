from __future__ import annotations
import importlib.util,json,os,subprocess,sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[2]
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,file);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module);return module
projection=load('bayestyper_projection',ROOT/'plugins/bayestyper/projection.py')
# Prevent another plugin's generic module name from changing this adapter's import.
sys.modules['projection']=projection
runner=load('bayestyper_runner',ROOT/'plugins/bayestyper/run.py')

def fixture(tmp_path, native_rows, panel_rows=None, sample='HG002'):
    ref=tmp_path/'ref.fa';ref.write_text('>chr1\nAAAAAAAAAA\n')
    Path(str(ref)+'.fai').write_text('chr1\t10\t6\t10\t11\n')
    panel=tmp_path/'panel.vcf'
    panel.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+(panel_rows or 'chr1\t3\tC1\tA\tAA\t99\tPASS\tPANGENOME_ALLELE_ID=C1\n'))
    native=tmp_path/'native.vcf';native.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n'+native_rows)
    return ref,panel,native,tmp_path/'calls.vcf'

@pytest.mark.parametrize('gt,expected',[('0/0','0/0'),('0/1','0/1'),('1/1','1/1'),('./1','./.')])
def test_preserves_native_quality_and_filter(tmp_path,gt,expected):
    ref,panel,native,out=fixture(tmp_path,f'chr1\t3\tC1\tA\tAA\t7.5\tAN0\t.\tGT:GQ\t{gt}:3\n')
    projection.project(native,panel,out,'HG002',reference=ref)
    fields=next(l for l in out.read_text().splitlines() if not l.startswith('#')).split('\t')
    assert fields[5:7]==['7.5','AN0']
    assert fields[9]==expected+':3'
    audit=json.loads(out.with_suffix('.projection.json').read_text())
    assert audit['quality_and_filter_preserved'] and audit['truth_used'] is False
    assert audit['reference_sha256'] == projection.digest(ref)

def test_multiallelic_dosage(tmp_path):
    ref,panel,native,out=fixture(tmp_path,'chr1\t3\tC1\tA\tAC,AA\t30\tPASS\t.\tGT\t2/2\n')
    projection.project(native,panel,out,'HG002',reference=ref)
    assert out.read_text().rstrip().endswith('GT\t1/1')

def test_missing_is_not_reference_and_wrong_id_is_not_evidence(tmp_path):
    ref,panel,native,out=fixture(tmp_path,'chr1\t3\tC1\tA\tAC\t30\tPASS\t.\tGT\t0/1\n')
    assert projection.project(native,panel,out,'HG002',reference=ref)[:2]==(0,1)
    assert 'PGBENCH_NATIVE_STATUS=not_observed' in out.read_text()
    assert out.read_text().rstrip().endswith('GT\t./.')

def test_ambiguous_mapping_never_chooses_highest_quality(tmp_path):
    ref,panel,native,out=fixture(tmp_path,'chr1\t3\ta\tA\tAA\t1\tPASS\t.\tGT\t0/1\nchr1\t3\tb\tA\tAA\t99\tPASS\t.\tGT\t1/1\n')
    projection.project(native,panel,out,'HG002',reference=ref)
    assert 'PGBENCH_NATIVE_STATUS=ambiguous_mapping' in out.read_text()
    assert out.read_text().rstrip().endswith('GT\t./.')

def test_rejects_wrong_native_sample(tmp_path):
    ref,panel,native,out=fixture(tmp_path,'',sample='OTHER')
    with pytest.raises(ValueError,match='single sample'):projection.project(native,panel,out,'HG002',reference=ref)

def test_kmc_enforces_budget_and_bayestyper_byte_counts():
    command=runner.kmc_command(Path('kmc'),Path('reads'),Path('prefix'),Path('tmp'),4,8192)
    assert '-t4' in command and '-m7' in command and '-sm' in command and '-cs255' in command
    assert not any(arg.startswith('-cx') for arg in command)
    with pytest.raises(ValueError):runner.kmc_command(Path('kmc'),Path('reads'),Path('prefix'),Path('tmp'),4,2048)


@pytest.mark.real_tools
def test_real_kmc_saturates_high_counts_and_retains_singletons(tmp_path):
    kmc, dump = os.environ.get('PGBENCH_REAL_KMC'), os.environ.get('PGBENCH_REAL_KMC_DUMP')
    if not kmc or not dump:
        pytest.skip('Set PGBENCH_REAL_KMC and PGBENCH_REAL_KMC_DUMP for real KMC regression')
    high = ('ACGTTGCA' * 7)[:55]
    singleton = ('GATCCAGT' * 7)[:55]
    fastq = tmp_path/'reads.fq'
    fastq.write_text(''.join(f'@r{i}\n{s}\n+\n'+('I'*55)+'\n'
                             for i,s in enumerate([high]*300+[singleton])))
    reads = tmp_path/'reads.list'; reads.write_text(str(fastq)+'\n')
    temporary = tmp_path/'tmp'; temporary.mkdir()
    observed = {}
    for recipe in ('legacy', 'fixed'):
        prefix = tmp_path/recipe
        command = runner.kmc_command(Path(kmc),reads,prefix,temporary,2,3072)
        if recipe == 'legacy':
            command[command.index('-cs255')] = '-cs65535'
        subprocess.run(command,check=True,capture_output=True)
        counts = tmp_path/(recipe+'.txt')
        subprocess.run([dump,str(prefix),str(counts)],check=True,capture_output=True)
        observed[recipe] = dict(line.split() for line in counts.read_text().splitlines())
    assert observed['legacy'].keys() == observed['fixed'].keys()
    assert sorted(map(int,observed['legacy'].values())) == [1,300]
    assert sorted(map(int,observed['fixed'].values())) == [1,255]

@pytest.mark.parametrize('sample',['../HG002','HG002\tM','x/y'])
def test_rejects_sample_path_injection(sample):
    with pytest.raises(ValueError):runner.validate_sample(sample)


@pytest.mark.parametrize('gt',['0/2','1/2','2/2'])
def test_spanning_deletion_is_not_a_complete_diploid_call(tmp_path,gt):
    ref,panel,native,out=fixture(tmp_path,f'chr1\t3\tC1\tA\tAA,*\t30\tPASS\t.\tGT\t{gt}\n')
    projection.project(native,panel,out,'HG002',reference=ref)
    assert out.read_text().rstrip().endswith('GT\t./.')


@pytest.mark.parametrize('contents',['', '##fileformat=VCFv4.2\n', 'chr1\t3\tC1\tA\tAA\t30\tPASS\t.\tGT\t0/1\n'])
def test_missing_native_header_cannot_fabricate_all_sites_no_calls(tmp_path,contents):
    ref,panel,native,out=fixture(tmp_path,'')
    native.write_text(contents)
    with pytest.raises(ValueError,match='header'):
        projection.project(native,panel,out,'HG002',reference=ref)
    assert not out.exists()


def test_reference_index_must_be_declared_and_adjacent(tmp_path):
    reference=tmp_path/'ref.fa';reference.write_text('>chr1\nA\n')
    index=tmp_path/'elsewhere.fai';index.write_text('chr1\t1\t6\t1\t2\n')
    with pytest.raises(ValueError,match='adjacent'):runner.validate_reference_index(reference,index)
    adjacent=Path(str(reference)+'.fai');adjacent.write_text(index.read_text())
    runner.validate_reference_index(reference,adjacent)


def test_same_position_candidates_merge_and_project_without_losing_ids(tmp_path):
    import gzip
    ref,panel,native,out=fixture(tmp_path,'',
        'chr1\t3\tC1\tA\tAA\t.\t.\tPANGENOME_ALLELE_ID=C1\n'
        'chr1\t3\tC2\tAAA\tA\t.\t.\tPANGENOME_ALLELE_ID=C2\n')
    prepared=tmp_path/'prepared.vcf.gz'
    stats=runner.prepare_candidates(panel,prepared)
    with gzip.open(prepared,'rt') as stream:
        row=next(line for line in stream if not line.startswith('#')).rstrip().split('\t')
    assert row[2:5]==['C1;C2','AAA','AAAA,A']
    assert stats==dict(candidate_count=2,native_site_count=1,merged_sites=1)
    native.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n'+
                     '\t'.join(row[:5]+['30','PASS','.','GT','1/2'])+'\n')
    assert projection.project(native,panel,out,'HG002',reference=ref)[:2]==(2,0)
    rows=[line.split('\t') for line in out.read_text().splitlines() if not line.startswith('#')]
    assert [(row[2],row[9]) for row in rows]==[('C1','0/1'),('C2','0/1')]


def test_unsorted_candidates_rejected_before_kmer_counting(tmp_path):
    ref,panel,native,out=fixture(tmp_path,'',
        'chr1\t4\tC1\tA\tAA\t.\t.\t.\nchr1\t3\tC2\tA\tAAA\t.\t.\t.\n')
    with pytest.raises(ValueError,match='sorted'):
        runner.prepare_candidates(panel,tmp_path/'prepared.vcf.gz')
