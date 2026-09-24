# PG-F1 is the sole scoring and release chain.
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
            units=PGF1_UNITS,
        output:
            query=PGF1_ROOT + "/{tool}/evaluation/evaluator-query.vcf.gz",
            audit=PGF1_ROOT + "/{tool}/evaluation/evaluator-query.audit.json",
        conda: "../envs/core.yaml"
        shell:
            "python workflow/scripts/materialize_pgf1_evaluator_query.py --linked-vcf {input.linked:q} --canonical-units {input.units:q} --output-vcf {output.query:q} --audit-json {output.audit:q}"

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

    rule finalize_pgf1_release:
        input:
            score=PGF1_ROOT + "/{tool}/evaluation/consensus/score.json",
            score_audit=PGF1_ROOT + "/{tool}/evaluation/consensus/pgf1_audit.json",
            ledger=PGF1_ROOT + "/{tool}/evaluation/consensus/pgf1_evidence.tsv.gz",
            canonical_units=PGF1_UNITS,
            unit_sha=PGF1_UNIT_SHA,
            all_sites=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/{{tool}}/canonical/all-sites.vcf.gz",
            query=PGF1_ROOT + "/{tool}/evaluation/evaluator-query.vcf.gz",
            query_audit=PGF1_ROOT + "/{tool}/evaluation/evaluator-query.audit.json",
            truvari_evidence=PGF1_ROOT + "/{tool}/evaluation/normalized/truvari.evidence.tsv.gz",
            aardvark_evidence=PGF1_ROOT + "/{tool}/evaluation/normalized/aardvark.evidence.tsv.gz",
            vcfdist_evidence=PGF1_ROOT + "/{tool}/evaluation/normalized/vcfdist.evidence.tsv.gz",
            truvari_raw=RESULTS_ROOT + "/evaluation/{tool}/truvari/votes.tsv",
            aardvark_raw=RESULTS_ROOT + "/evaluation/{tool}/aardvark/votes.tsv",
            vcfdist_raw=RESULTS_ROOT + "/evaluation/{tool}/vcfdist/votes.tsv",
            invariants=PGF1_ROOT + "/leaderboard-invariants.json",
            admission=PGF1_ROOT + "/{tool}/evaluation/leaderboard-admission.json",
            context=RESULTS_ROOT + "/provenance/run-context.json",
            config=CONFIG_PATH,
            rule_source="workflow/rules/pgf1.smk",
            rule_executor=RULE_EXECUTOR,
            script="workflow/scripts/finalize_pgf1_release.py",
            environment=CORE_ENV_SPEC,
            provenance_library=PROVENANCE_LIBRARY,
            link_manifest=lambda wc: semantic_rule_manifest(
                "link_pangenome_alleles", f"{SAMPLE_ID}.{wc.tool}.{OFFICIAL_MODE}"
            ),
            truvari_manifest=RESULTS_ROOT + "/provenance/rules/evaluate_truvari/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
            aardvark_manifest=RESULTS_ROOT + "/provenance/rules/evaluate_aardvark/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
            vcfdist_manifest=RESULTS_ROOT + "/provenance/rules/evaluate_vcfdist/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
        output:
            package=PGF1_ROOT + "/{tool}/release/score-package.json",
            audit=PGF1_ROOT + "/{tool}/release/provenance-audit.json",
            seal=PGF1_ROOT + "/{tool}/release/seal.json",
            rule_manifest=(
                RESULTS_ROOT + "/provenance/rules/finalize_pgf1_release/"
                f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
            ),
        log:
            LOG_ROOT + "/rules/finalize_pgf1_release/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".log",
        benchmark:
            BENCHMARK_ROOT + "/rules/finalize_pgf1_release/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".jsonl",
        conda:
            "../envs/core.yaml"
        shell:
            """
            {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \\
              --rule-name finalize_pgf1_release \\
              --job-key {SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE} \\
              --run-id {RUN_ID:q} --module-or-tool-id pgbench-core \\
              --snakefile-path {input.rule_source:q} --rule-source-path {input.rule_source:q} \\
              --script-or-wrapper-path {input.script:q} --config-snapshot {input.config:q} \\
              --score-profile {config[catalogs][score_weights]:q} --run-context {input.context:q} \\
              --reference {config[reference][fasta]:q} --truth-profile {config[truth][primary]:q} \\
              --snakemake-version {SNAKEMAKE_VERSION:q} --execution-profile local \\
              --random-seed {config[execution][random_seed]} --threads 1 --resource mem_mb=2048 \\
              --wildcard tool={wildcards.tool:q} \\
              --input {input.score:q} --input {input.score_audit:q} --input {input.ledger:q} \\
              --input {input.canonical_units:q} --input {input.unit_sha:q} --input {input.all_sites:q} \\
              --input {input.query:q} --input {input.query_audit:q} \\
              --input {input.truvari_evidence:q} --input {input.aardvark_evidence:q} --input {input.vcfdist_evidence:q} \\
              --input {input.truvari_raw:q} --input {input.aardvark_raw:q} --input {input.vcfdist_raw:q} \\
              --input {input.invariants:q} --input {input.admission:q} --input {input.context:q} \\
              --input {input.config:q} --input {input.rule_source:q} --input {input.rule_executor:q} \\
              --input {input.script:q} --input {input.environment:q} --input {input.provenance_library:q} \\
              --output {output.package:q} --output {output.audit:q} --output {output.seal:q} \\
              --upstream-manifest {input.link_manifest:q} --upstream-manifest {input.truvari_manifest:q} \\
              --upstream-manifest {input.aardvark_manifest:q} --upstream-manifest {input.vcfdist_manifest:q} \\
              --manifest-output {output.rule_manifest:q} -- \\
              {PYTHON_EXECUTABLE:q} {input.script:q} \\
                --score {input.score:q} --score-audit {input.score_audit:q} --ledger {input.ledger:q} \\
                --canonical-units {input.canonical_units:q} --unit-sha-file {input.unit_sha:q} \\
                --all-sites {input.all_sites:q} --query {input.query:q} --query-audit {input.query_audit:q} \\
                --truvari-evidence {input.truvari_evidence:q} --aardvark-evidence {input.aardvark_evidence:q} --vcfdist-evidence {input.vcfdist_evidence:q} \\
                --raw-truvari-ledger {input.truvari_raw:q} --raw-aardvark-ledger {input.aardvark_raw:q} --raw-vcfdist-ledger {input.vcfdist_raw:q} \\
                --invariants {input.invariants:q} --admission {input.admission:q} \\
                --run-context {input.context:q} --config {input.config:q} \\
                --output-package {output.package:q} --output-audit {output.audit:q} --output-seal {output.seal:q} \\
              > {log:q} 2>&1
            """
