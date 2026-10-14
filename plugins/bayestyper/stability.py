"""Truth-blind noise screening with explicit successful-attempt binding."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import statistics

EXPECTED_CHAINS, BURN, DRAWS = 20, 100, 250
RHAT_FLAG = 1.05  # Local screening flag, not an upstream requirement.

def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def identity(p):
    s = p.stat()
    return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)


def split_rhat(chains):
    length = len(chains[0]) // 2
    splits = [half for chain in chains for half in (chain[:length], chain[-length:])]
    means = [statistics.mean(x) for x in splits]
    within = statistics.mean(statistics.variance(x) for x in splits)
    between = length * statistics.variance(means)
    if within == 0:
        return None, 'DEGENERATE_CHAINS_IDENTICAL' if between == 0 else 'DEGENERATE_CHAINS_DISAGREE'
    return math.sqrt(((length - 1) * within / length + between / length) / within), None


def analyze_file(path, *, noise_genotyping=False):
    before = identity(path)
    text = path.read_text()
    require(identity(path) == before, 'Noise file is changing; check again after it is complete')
    lines = text.splitlines()
    if not lines:
        return {'status': 'PENDING', 'reason': 'Empty noise output'}
    header = lines[0].split('\t')
    require(header[:2] == ['Chain', 'Iteration'] and len(header) >= 3, 'Unexpected noise table header')
    samples = header[2:]
    require(len(set(samples)) == len(samples), 'Duplicate sample columns')
    rows, final = {}, None
    for line in lines[1:]:
        f = line.split('\t')
        # The writer may have only partially written the last row.
        if len(f) != len(header):
            return {'status': 'PENDING', 'reason': 'Incomplete row; do not assess a file still being written'}
        chain, iteration = int(f[0]), int(f[1])
        values = list(map(float, f[2:]))
        require(all(math.isfinite(v) and v >= 0 for v in values), 'Nonfinite/negative noise rates')
        if chain == 0:
            require(iteration == 0 and final is None, 'Duplicate/unexpected final summary')
            final = values
            continue
        require(1 <= chain <= EXPECTED_CHAINS and 0 <= iteration <= BURN + DRAWS, 'Unexpected chain/iteration')
        require((chain, iteration) not in rows, 'Duplicate chain iteration')
        rows[chain, iteration] = values
    expected = {(c, i) for c in range(1, EXPECTED_CHAINS + 1) for i in range(0, BURN + DRAWS + 1)}
    # v1.5 joint mode has no sequential-mode final-average 0/0 row.
    if set(rows) != expected or (final is None and not noise_genotyping):
        return {'status': 'PENDING', 'rows_present': len(rows), 'expected_rows': len(expected),
                'final_summary_present': final is not None, 'reason': 'Waiting for all 20 complete sampling chains'}
    summary = {}
    for column, sample in enumerate(samples):
        chains = [[rows[c, i][column] for i in range(BURN + 1, BURN + DRAWS + 1)]
                  for c in range(1, EXPECTED_CHAINS + 1)]
        means = [statistics.mean(c) for c in chains]
        pooled_mean = statistics.mean(means)
        rhat, degenerate = split_rhat(chains)
        summary[sample] = {
            'retained_chains': EXPECTED_CHAINS, 'burn_in_discarded_per_chain': BURN,
            'draws_per_chain': DRAWS, 'chain_means': means, 'mean_noise_rate': pooled_mean,
            'chain_mean_min': min(means), 'chain_mean_max': max(means),
            'chain_mean_cv': statistics.stdev(means) / pooled_mean if pooled_mean else None,
            'classical_split_rhat': rhat, 'degenerate_chains': degenerate,
            'native_final_mean': final[column] if final is not None else None,
            # Native table uses limited decimal precision, so do not demand
            # exact equality to its independently accumulated final average.
            'final_mean_absolute_difference': abs(final[column] - pooled_mean) if final is not None else None,
            'first_half_mean': statistics.mean(v for c in chains for v in c[:DRAWS // 2]),
            'second_half_mean': statistics.mean(v for c in chains for v in c[-(DRAWS // 2):]),
            'screen_flag': degenerate is not None or rhat is None or rhat > RHAT_FLAG,
        }
    return {'status': 'REVIEW_REQUIRED' if any(r['screen_flag'] for r in summary.values()) else 'NO_SCREEN_FLAG',
            'path': str(path), 'sha256': sha(path), 'samples': summary,
            'interpretation': 'Screening only. No-screen-flag does not prove mixing, genotype accuracy or unbiased noise estimation.'}


def inspect_run(acceptance, attempt_id):
    """Bind a successful attempt explicitly, never select using its score."""
    base = Path(acceptance).resolve()
    require(attempt_id and all(c in '0123456789abcdef' for c in attempt_id), 'Invalid attempt ID')
    spec = json.loads((base / 'run-spec.json').read_text())
    config = Path(spec['config'])
    require(sha(config) == spec['configuration_sha256'], 'Frozen configuration changed')
    root = Path(spec['result']) / 'HG002/end_to_end_from_reads/bayestyper'
    require(root.resolve().is_relative_to(base), 'Result is outside acceptance directory')
    attempt = json.loads((root / 'meta' / f'attempt.{attempt_id}.json').read_text())
    require(attempt['attempt_id'] == attempt_id and attempt['run_id'] == spec['run_id'], 'Attempt identity mismatch')
    require(attempt['status'] == 'success' and attempt['exit_code'] == 0, 'Select a completed successful attempt')
    require(attempt['sample_id'] == 'HG002', 'Unexpected sample')
    output = root / 'raw/calls.vcf'
    require(sha(output) == attempt['output']['sha256'], 'Published output differs from selected attempt')
    audit = json.loads(output.with_suffix('.native-audit.json').read_text())
    require(audit['output_sha256'] == sha(output), 'Native audit output binding mismatch')
    work = root / '.pgbench_attempts' / attempt_id / 'output/bayestyper-work'
    native = Path(audit['artifacts']['native_vcf']['path'])
    require(native.resolve() == (work / 'native.vcf.gz').resolve(), 'Native VCF belongs to another attempt')
    require(sha(native) == audit['artifacts']['native_vcf']['sha256'], 'Native VCF changed')
    frozen = json.loads((base / 'manifests/source-state.json').read_text())['implementation_hashes']
    source_binding = []
    # Today's unrelated tools need not match this historical execution.
    # BayesTyper files must match their executed hashes, in place or archived.
    for source, expected in frozen.items():
        if '/plugins/bayestyper/' not in source:
            continue
        relative = source.split('/plugins/bayestyper/', 1)[1]
        alternatives = [Path(source), base / 'manifests/adapter-source-snapshot' / relative]
        alternatives.extend(base.glob('*/executed-source-before-metadata-fix/' + relative))
        alternatives.extend(base.glob('*/executed-source-before-stability-fix/' + relative))
        bound = next((p for p in alternatives if p.is_file() and sha(p) == expected), None)
        require(bound is not None, 'Executed source not preserved: ' + source)
        source_binding.append({'source': source, 'executed_sha256': expected, 'verified_at': str(bound)})
    require(source_binding, 'No executed BayesTyper source bindings')
    require(frozen.get(attempt['runner_path']) == attempt['runner_sha256'], 'Attempt runner differs from frozen source')
    units = sorted(work.glob('bayestyper_unit_*/variant_clusters.bin'))
    require(units, 'Successful attempt has no genotyping units')
    joint = audit.get('parameter_policy', {}).get('noise_genotyping', False)
    require(isinstance(joint, bool), 'Invalid noise_genotyping policy')
    diagnostics = []
    noise_bound_at_execution = True
    for unit in units:
        noise = unit.parent / 'calls_noise_parameters.txt'
        require(noise.is_file() and not noise.is_symlink(), 'Missing/redirected noise output: ' + str(noise))
        row = analyze_file(noise, noise_genotyping=joint)
        require(set(row.get('samples', {})) == {'HG002'}, 'Incomplete noise output or unexpected sample')
        recorded = audit['artifacts'].get('noise_parameters_' + unit.parent.name)
        if recorded:
            require(Path(recorded['path']).resolve() == noise.resolve() and recorded['sha256'] == row['sha256'], 'Noise hash binding mismatch')
        else:
            noise_bound_at_execution = False
        diagnostics.append(row)
    return {'schema': 'pgbench.external-bayestyper-noise-screen.v2', 'run_id': spec['run_id'],
            'attempt_id': attempt_id, 'created_utc': datetime.now(timezone.utc).isoformat(),
            'status': 'REVIEW_REQUIRED' if any(d['status'] == 'REVIEW_REQUIRED' for d in diagnostics) else 'NO_SCREEN_FLAG',
            'noise_genotyping': joint, 'noise_bound_at_execution': noise_bound_at_execution,
            'unit_diagnostics': diagnostics, 'executed_source_bindings': source_binding,
            'configuration_sha256': sha(config), 'native_sha256': sha(native), 'output_sha256': sha(output),
            'program_sha256': sha(__file__), 'writes_to_caller_or_score': False,
            'method': {'statistic': 'classical split R-hat', 'local_review_threshold': RHAT_FLAG,
                       'threshold_is_upstream_requirement': False, 'chains': EXPECTED_CHAINS, 'burn_in': BURN, 'draws': DRAWS},
            'limitations': ['No rank-normalized R-hat or ESS; no convergence or genotype accuracy proof.',
                            'Historical unbound noise hashes cannot be retroactively authenticated at execution time.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--acceptance', required=True, type=Path)
    parser.add_argument('--attempt', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    destination = args.output.resolve()
    allowed = args.acceptance.resolve() / 'supplemental-noise-stability'
    require(destination.parent == allowed and not destination.exists(), 'Use a new immediate child of supplemental-noise-stability')
    result = inspect_run(args.acceptance, args.attempt)
    destination.mkdir(parents=True, exist_ok=False)
    (destination / 'diagnostics.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (destination / 'report.md').write_text('# BayesTyper noise screening\n\nStatus: **' + result['status'] + '**\n\n'
        + 'Run: `' + result['run_id'] + '`; attempt: `' + result['attempt_id'] + '`.\n\n'
        + 'Diagnostic only; caller outputs, score and seal are unchanged. No-screen-flag is not proof of convergence.\n\n'
        + 'See diagnostics.json for source binding, file hashes, chain means and classical split R-hat.\n', encoding='utf-8')
    print(json.dumps({'status': result['status'], 'report': str(destination / 'report.md'), 'run_modified': False}))


if __name__ == '__main__':
    main()


