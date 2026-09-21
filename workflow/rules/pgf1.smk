# PG-F1 is a separate scoring/release chain.  It is intentionally selected only
# by the pgbench_pgf1_v1 profile, so frozen ME-F1 releases retain their exact
# historical DAG and provenance.
if PGF1_ENABLED:
    PGF1_ROOT = f"{RESULTS_ROOT}/{SAMPLE_ID}/{config['benchmark_contract']['track']}"
    PGF1_UNITS = f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/{SAMPLE_ID}.canonical_units.tsv.gz"
    PGF1_UNIT_SHA = f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/{SAMPLE_ID}.canonical_unit_set.sha256"
    PGF1_UNIT_AUDIT = f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/{SAMPLE_ID}.canonical_units.audit.json"

    rule materialize_canonical_units:
        input:
            panel=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/{SAMPLE_ID}.blinded.vcf",
            hidden=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/{SAMPLE_ID}.hidden.tsv",
        output:
            units=PGF1_UNITS, sha=PGF1_UNIT_SHA, audit=PGF1_UNIT_AUDIT,
        conda: "../envs/core.yaml"
        shell:
            "python workflow/scripts/materialize_canonical_units.py --canonical-panel {input.panel:q} --hidden-truth-ledger {input.hidden:q} --output-units {output.units:q} --output-sha256 {output.sha:q} --audit-json {output.audit:q}"

    rule materialize_pgf1_evaluator_query:
        input:
            linked=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/{{tool}}/canonical/linked.vcf",
        output:
            query=PGF1_ROOT + "/{tool}/evaluation/evaluator-query.vcf.gz",
        conda: "../envs/core.yaml"
        shell:
            "python workflow/scripts/materialize_pgf1_evaluator_query.py --linked-vcf {input.linked:q} --output-vcf {output.query:q}"

    rule normalize_pgf1_evidence:
        input:
            units=PGF1_UNITS,
            all_sites=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/{{tool}}/canonical/all-sites.vcf.gz",
            ledger=RESULTS_ROOT + "/evaluation/{tool}/{evaluator}/votes.tsv",
        output:
            evidence=PGF1_ROOT + "/{tool}/evaluation/normalized/{evaluator}.evidence.tsv.gz",
        params:
            release=lambda wc: f"{RUN_ID}-PGF1",
        conda: "../envs/core.yaml"
        wildcard_constraints:
            evaluator="truvari|aardvark|vcfdist"
        shell:
            "python workflow/scripts/normalize_evaluator_evidence.py --units {input.units:q} --all-sites {input.all_sites:q} --raw-ledger {input.ledger:q} --evaluator {wildcards.evaluator:q} --release-id {params.release:q} --run-id {RUN_ID:q} --sample-id {SAMPLE_ID:q} --tool-id {wildcards.tool:q} --output {output.evidence:q}"

    rule materialize_pg_f1:
        input:
            units=PGF1_UNITS,
            all_sites=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/{{tool}}/canonical/all-sites.vcf.gz",
            truvari=PGF1_ROOT + "/{tool}/evaluation/normalized/truvari.evidence.tsv.gz",
            aardvark=PGF1_ROOT + "/{tool}/evaluation/normalized/aardvark.evidence.tsv.gz",
            vcfdist=PGF1_ROOT + "/{tool}/evaluation/normalized/vcfdist.evidence.tsv.gz",
        output:
            ledger=PGF1_ROOT + "/{tool}/evaluation/consensus/pgf1_evidence.tsv.gz",
            audit=PGF1_ROOT + "/{tool}/evaluation/consensus/pgf1_audit.json",
            score=PGF1_ROOT + "/{tool}/evaluation/consensus/score.json",
        conda: "../envs/core.yaml"
        shell:
            "python workflow/scripts/materialize_pgf1_evidence.py --units {input.units:q} --all-sites {input.all_sites:q} --truvari-evidence {input.truvari:q} --aardvark-evidence {input.aardvark:q} --vcfdist-evidence {input.vcfdist:q} --ledger-output {output.ledger:q} --audit-output {output.audit:q} --score-output {output.score:q} --run-id {RUN_ID:q} --sample-id {SAMPLE_ID:q} --tool-id {wildcards.tool:q} --track {config[benchmark_contract][track]:q}"

    rule materialize_pgf1_invariants:
        input:
            panel=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/{SAMPLE_ID}.blinded.vcf",
            units=PGF1_UNITS,
            config=CONFIG_PATH,
            evaluator=config["catalogs"]["evaluator_profile"],
        output:
            invariants=PGF1_ROOT + "/leaderboard-invariants.json",
        conda: "../envs/core.yaml"
        shell:
            "python workflow/scripts/materialize_pgf1_invariants.py --config {input.config:q} --panel {input.panel:q} --units {input.units:q} --evaluator-profile {input.evaluator:q} --output {output.invariants:q}"

    rule validate_pgf1_leaderboard_admission:
        input:
            score=PGF1_ROOT + "/{tool}/evaluation/consensus/score.json",
            invariants=PGF1_ROOT + "/leaderboard-invariants.json",
        output:
            admission=PGF1_ROOT + "/{tool}/evaluation/leaderboard-admission.json",
        conda: "../envs/core.yaml"
        shell:
            "python workflow/scripts/validate_leaderboard_admission.py --score {input.score:q} --invariants {input.invariants:q} --output {output.admission:q}"
