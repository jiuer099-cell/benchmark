import gzip
import json
import pathlib
import tempfile
import unittest
import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / 'plugins/gatk-sv'))
import run_rg_compatible as adapter
import resume


class CompressedEvidenceTests(unittest.TestCase):
    def fixture(self, folder, counts, compressed=True):
        outputs = {}
        for name in ['pesr_disc', 'pesr_split', 'pesr_sd']:
            path = folder / (name + '.txt.gz')
            with gzip.open(path, 'wt') as stream:
                stream.write('chr1\t1\t2\tHG002\n')
            outputs['GatherSampleEvidence.' + name] = str(path)
        path = folder / ('counts.tsv.gz' if compressed else 'counts.tsv')
        with (gzip.open(path, 'wt') if compressed else path.open('w')) as stream:
            stream.write(counts)
        outputs['GatherSampleEvidence.coverage_counts'] = str(path)
        return dict(outputs=outputs)

    def test_actual_gzip_counts_and_plain_counts_both_pass(self):
        for compressed in [True, False]:
            with self.subTest(compressed=compressed), tempfile.TemporaryDirectory() as d:
                m = self.fixture(pathlib.Path(d), '@HD\tVN:1.6\nCONTIG\tSTART\tEND\tCOUNT\nchr1\t1\t100\t42\n', compressed)
                self.assertEqual(adapter.evidence_gate(m)['status'], 'PASS')

    def test_all_zero_counts_still_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            m = self.fixture(pathlib.Path(d), 'CONTIG\tSTART\tEND\tCOUNT\nchr1\t1\t100\t0\n')
            with self.assertRaisesRegex(ValueError, 'All-zero'):
                adapter.evidence_gate(m)

    def test_corrupt_gzip_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            m = self.fixture(pathlib.Path(d), 'chr1\t1\t100\t42\n')
            pathlib.Path(m['outputs']['GatherSampleEvidence.coverage_counts']).write_bytes(b'not gzip data')
            with self.assertRaises(gzip.BadGzipFile):
                adapter.evidence_gate(m)

    def test_altered_checkpoint_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = pathlib.Path(d)
            historical = root / 'old'
            historical.mkdir()
            meta = historical / 'metadata.json'
            meta.write_text('original')
            manifest = root / 'resume-inputs.json'
            manifest.write_text(json.dumps(dict(resume_from=str(historical), sha256={str(meta):adapter.base.digest(meta)})))
            config = dict(output_dir=str(root / 'new'), resume_manifest=str(manifest), resume_readonly_paths=[str(historical)])
            meta.write_text('changed')
            with self.assertRaisesRegex(ValueError, 'checkpoint changed'):
                resume.authenticate(config)


if __name__ == '__main__':
    unittest.main()
