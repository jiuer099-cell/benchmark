#!/usr/bin/env python3
"""Resume official GATK-SV using authenticated completed sample evidence."""
import argparse
import json
import os
import pathlib
import subprocess
import datetime
import run_rg_compatible as adapter

ROOT = pathlib.Path(__file__).resolve().parent
PREFIX = 'GATKSVPipelineSingleSample.'
EVIDENCE_MAP = {
    'case_manta_vcf': 'manta_vcf', 'case_wham_vcf': 'wham_vcf',
    'case_scramble_vcf': 'scramble_vcf', 'case_counts_file': 'coverage_counts',
    'case_pe_file': 'pesr_disc', 'case_sr_file': 'pesr_split', 'case_sd_file': 'pesr_sd',
}


def authenticate(config):
    manifest_path = pathlib.Path(config['resume_manifest'])
    manifest = json.loads(manifest_path.read_text())
    old = pathlib.Path(manifest['resume_from']).resolve()
    out = pathlib.Path(config['output_dir']).resolve()
    if out == old or out.is_relative_to(old) or old.is_relative_to(out):
        raise ValueError('Recovery output overlaps historical attempt')
    if config['resume_readonly_paths'] != [str(old)]:
        raise ValueError('Historical attempt must be mounted read-only')
    for path, expected in manifest['sha256'].items():
        if adapter.base.digest(path) != expected:
            raise ValueError('Recovery checkpoint changed: ' + path)
    previous = json.loads((old / 'attempt.json').read_text())
    new_settings = {k: v for k, v in config.items() if k not in {'output_dir', 'resume_manifest', 'resume_readonly_paths'}}
    old_settings = {k: v for k, v in previous['config'].items() if k != 'output_dir'}
    if old_settings != new_settings:
        raise ValueError('Recovery settings differ from successful evidence production')
    if previous.get('status') != 'failed' or "'utf-8' codec can't decode byte 0x8b" not in previous.get('error', ''):
        raise ValueError('Not the authenticated gzip-validation failure checkpoint')
    for key, identity in [('bam_or_cram', 'source_identity'), ('alignment_index', 'source_index_identity')]:
        if adapter.rg_compat.identity(pathlib.Path(config[key])) != previous[identity]:
            raise ValueError('Frozen source alignment changed: ' + key)
    audit = previous['rg_compatibility_audit']
    if audit['status'] != 'PASS' or not audit['every_read_has_expected_rg'] or not audit['every_alignment_field_and_non_RG_aux_preserved']:
        raise ValueError('Prior read-group compatibility audit failed')
    if audit['normalized_records_sha256'] != config['source_alignment_records_sha256']:
        raise ValueError('Read-group conversion differs from frozen source records')
    alignment = pathlib.Path(audit['output']).resolve()
    index = pathlib.Path(str(alignment) + '.bai')
    for p in [alignment, index]:
        p.relative_to(old)
        if adapter.rg_compat.identity(p) != manifest['alignment_identity'][str(p)]:
            raise ValueError('Historical compatible alignment changed: ' + str(p))
    metadata = {}
    for part in ['manta', 'evidence']:
        m = json.loads((old / part / 'metadata.json').read_text())
        if m['status'] != 'Succeeded':
            raise ValueError('Cannot resume unfinished ' + part)
        stage_plan = json.loads((old / part / 'stage-plan.json').read_text())
        if adapter.base.digest(stage_plan['workflow']) != stage_plan['workflow_sha256']:
            raise ValueError('Native WDL changed: ' + part)
        if adapter.base.digest(old / part / 'inputs.json') != stage_plan['inputs_sha256']:
            raise ValueError('Native stage inputs changed: ' + part)
        metadata[part] = m
    for key in EVIDENCE_MAP.values():
        path = pathlib.Path(metadata['evidence']['outputs']['GatherSampleEvidence.' + key]).resolve()
        path.relative_to(old)
        if str(path) not in manifest['sha256']:
            raise ValueError('Unauthenticated native evidence: ' + key)
    return old, previous, alignment, index, metadata


def prepare(config):
    out = adapter.validate_inputs(config)
    old, previous, alignment, index, metadata = authenticate(config)
    original = pathlib.Path(config['bam_or_cram'])
    header = subprocess.check_output([config['samtools'], 'view', '-H', str(original)], text=True)
    adapter.rg_compat.single_rg(header, config['sample_id'])
    ref_check = adapter.rg_compat.original_reference_check(
        header, pathlib.Path(config['alignment_reference_dict']).read_text(), adapter.base.parse_dictionary)
    if adapter.base.digest(config['alignment_freeze_manifest']) != config['alignment_freeze_manifest_sha256']:
        raise ValueError('Frozen source manifest changed')
    for bam, bai in [(original, pathlib.Path(config['alignment_index'])), (alignment, index)]:
        subprocess.run([config['samtools'], 'quickcheck', '-v', str(bam)], check=True)
        subprocess.run([config['samtools'], 'view', '-c', '-X', str(bam), str(bai), 'chr1:1-1'], check=True, stdout=subprocess.DEVNULL)
    primary, workflow = adapter.base.render_inputs(pathlib.Path(config['download_root']), config['sample_id'], alignment, index)
    adapter.base.compatible_header(header, pathlib.Path(primary[PREFIX + 'reference_dict']).read_text(), config['sample_id'],
                                  alignment_reference_dict=pathlib.Path(config['alignment_reference_dict']).read_text())
    schema = adapter.input_schema(config, workflow.parent / 'GatherSampleEvidence.wdl')
    expected_gse = adapter.map_gse_inputs(primary, schema, metadata['manta']['outputs'], alignment, index)
    if expected_gse != json.loads((old / 'evidence' / 'inputs.json').read_text()):
        raise ValueError('Official evidence configuration no longer matches completed stage')
    proof = adapter.evidence_gate(metadata['evidence'])
    for target, key in EVIDENCE_MAP.items():
        primary[PREFIX + target] = metadata['evidence']['outputs']['GatherSampleEvidence.' + key]
    return out, previous, alignment, index, primary, workflow, proof, ref_check


def main(argv=None):
    if os.environ.get('PGBENCH_RESOLVED_INPUTS'):
        raise ValueError('Discovery is unsupported in fixed-panel PG-F1')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=pathlib.Path, required=True)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text())
    out = adapter.validate_inputs(config)
    plan = None
    if args.execute:
        out.mkdir(parents=True, exist_ok=False)
        plan = dict(config=config, status='running', resumed_from=json.loads(pathlib.Path(config['resume_manifest']).read_text())['resume_from'],
                    adapter_sha256=adapter.base.digest(__file__), native_filters='unchanged', truth_used=False,
                    primary_fixed_panel_score_eligible=False, skipped_stages=['RG_conversion', 'Manta', 'GatherSampleEvidence'])
        adapter.save(out / 'attempt.json', plan)
        adapter.status(out, 'validating_completed_evidence')
    try:
        out, previous, alignment, index, primary, workflow, proof, ref_check = prepare(config)
        if not args.execute:
            result = dict(status='PREFLIGHT_PASS', evidence_validation=proof, reference_check=ref_check,
                          native_filters='unchanged', skipped_stages=['RG_conversion', 'Manta', 'GatherSampleEvidence'])
            adapter.save(ROOT / 'preflight-report.json', result)
            print(json.dumps(result), flush=True)
            return 0
        adapter.save(out / 'evidence-validation.json', proof)
        plan['source_identity'] = previous['source_identity']
        plan['source_index_identity'] = previous['source_index_identity']
        plan['rg_compatibility_audit'] = previous['rg_compatibility_audit']
        adapter.save(out / 'attempt.json', plan)
        adapter.status(out, 'official_single_sample_downstream', reused_evidence=True)
        meta = adapter.run_stage(config, out / 'pipeline', 'GATKSVPipelineSingleSample', workflow, primary, alignment, index)
        plan['native_output'] = adapter.base.publish_native(meta, out / 'pipeline', config['sample_id'])
        plan['status'] = 'succeeded'
        for key, identity in [('bam_or_cram', 'source_identity'), ('alignment_index', 'source_index_identity')]:
            if adapter.rg_compat.identity(pathlib.Path(config[key])) != plan[identity]:
                raise ValueError('Frozen source alignment changed during downstream execution')
        adapter.status(out, 'complete', native_output=plan['native_output'])
    except BaseException as exc:
        if plan is not None:
            plan.update(status='failed', error=str(exc))
            adapter.status(out, 'failed', error=str(exc))
        raise
    finally:
        if plan is not None:
            plan['finished_at_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
            adapter.save(out / 'attempt.json', plan)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
