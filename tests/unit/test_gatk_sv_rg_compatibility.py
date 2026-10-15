import hashlib,importlib.util,json,pathlib,subprocess,sys,unittest
ROOT=pathlib.Path(__file__).resolve().parent/'gatk-sv'
if not ROOT.is_dir():ROOT=pathlib.Path(__file__).resolve().parents[2]/'plugins/gatk-sv'
sys.path.insert(0,str(ROOT))
import run_rg_compatible as runner
import rg_compat


class CompatibilityTests(unittest.TestCase):
    def test_single_declared_group(self):
        self.assertEqual(rg_compat.single_rg('@RG\tID:r\tSM:HG002','HG002'),'r')
        for h in ('','@RG\tID:r\tSM:OTHER','@RG\tID:r\tSM:HG002\n@RG\tID:s\tSM:HG002'):
            with self.assertRaises(ValueError):rg_compat.single_rg(h,'HG002')

    def test_manta_requires_exact_reference_not_superset(self):
        h='@SQ\tSN:chr1\tLN:10\tM5:a'
        self.assertEqual(rg_compat.original_reference_check(h,h,runner.base.parse_dictionary)['contigs'],1)
        with self.assertRaises(ValueError):rg_compat.original_reference_check(h,h+'\n@SQ\tSN:HLA\tLN:9',runner.base.parse_dictionary)

    def test_official_input_routing(self):
        primary={'GATKSVPipelineSingleSample.gatk_docker':'pinned','GATKSVPipelineSingleSample.scramble_docker':'scramble',
                 'GATKSVPipelineSingleSample.truth':'NOT_ALLOWED'}
        schema={'GatherSampleEvidence.gatk_docker':None,'GatherSampleEvidence.scramble_docker':None}
        out=runner.map_gse_inputs(primary,schema,{'Manta.vcf':'native.vcf.gz','Manta.index':'native.vcf.gz.tbi'},'compat.bam','compat.bai')
        self.assertEqual(out['GatherSampleEvidence.gatk_docker'],'pinned')
        self.assertEqual(out['GatherSampleEvidence.scramble_docker'],'scramble')
        self.assertIsNone(out['GatherSampleEvidence.manta_docker'])
        self.assertFalse(out['GatherSampleEvidence.run_localize_reads'])
        self.assertNotIn('GatherSampleEvidence.truth',out)

    def test_auditor_preserves_all_non_rg_fields(self):
        record=b'read1\t99\tchr1\t4\t60\t4M\t=\t9\t9\tACGT\tIIII\tNM:i:0\tAS:i:4\n'
        augmented=record[:-1]+b'\tRG:Z:r\n'
        result=subprocess.run([str(ROOT/'rg_record_audit'),'r'],input=augmented,stdout=subprocess.PIPE,check=True)
        audit=json.loads(result.stdout)
        self.assertEqual(audit['normalized_records_sha256'],hashlib.sha256(record).hexdigest())
        self.assertEqual(audit['records'],1)

    def test_auditor_rejects_missing_wrong_duplicate_rg(self):
        base=b'read1\t99\tchr1\t4\t60\t4M\t=\t9\t9\tACGT\tIIII'
        for suffix in (b'\n',b'\tRG:Z:other\n',b'\tRG:Z:r\tRG:Z:r\n'):
            result=subprocess.run([str(ROOT/'rg_record_audit'),'r'],input=base+suffix,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            self.assertNotEqual(result.returncode,0)

    def test_auditor_does_not_mask_non_rg_changes(self):
        a=b'read1\t99\tchr1\t4\t60\t4M\t=\t9\t9\tACGT\tIIII\tRG:Z:r\n'
        b=a.replace(b'\t60\t',b'\t59\t')
        hashes=[json.loads(subprocess.check_output([str(ROOT/'rg_record_audit'),'r'],input=x))['normalized_records_sha256'] for x in (a,b)]
        self.assertNotEqual(*hashes)

    def test_fixed_panel_invocation_refused(self):
        import os
        os.environ['PGBENCH_RESOLVED_INPUTS']='x'
        try:
            with self.assertRaises(ValueError):runner.main([])
        finally:os.environ.pop('PGBENCH_RESOLVED_INPUTS')


if __name__=='__main__':unittest.main()
