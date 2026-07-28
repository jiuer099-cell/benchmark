FORMAL_EVALUATORS = ("truvari", "aardvark", "vcfdist")

if not SYNTHETIC_MODE:
    rule run_formal_evaluator:
        input:
            query=(
                f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/"
                "{tool}/canonical/linked.vcf"
            ),
            link_rule_manifest=(
                "results/provenance/rules/link_pangenome_alleles/"
                f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
            ),
            truth=config["evaluation"]["truth_vcf"],
            truth_index=config["evaluation"]["truth_vcf"] + ".tbi",
            regions=config["evaluation"]["benchmark_bed"],
            reference=config["reference"]["fasta"],
            pangenome_manifest=(
                f"results/pangenome/{PANGENOME_ID}/manifest.yaml"
            ),
            context="results/provenance/run-context.json",
            config=CONFIG_PATH,
            score_profile=config["catalogs"]["score_weights"],
            rule_source="workflow/rules/evaluation.smk",
            rule_executor=RULE_EXECUTOR,
            script="workflow/scripts/run_formal_evaluator.py",
            environment=CORE_ENV_SPEC,
            provenance_library=PROVENANCE_LIBRARY,
        output:
            ledger="results/evaluation/{tool}/{evaluator}/votes.tsv",
            complete="results/evaluation/{tool}/{evaluator}/.complete.json",
            rule_manifest=(
                "results/provenance/rules/evaluate_{evaluator}/"
                f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
            ),
        log:
            (
                f"logs/rules/evaluate_formal/{SAMPLE_ID}."
                "{tool}.{evaluator}." + OFFICIAL_MODE + ".log"
            ),
        benchmark:
            (
                f"benchmarks/rules/evaluate_formal/{SAMPLE_ID}."
                "{tool}.{evaluator}." + OFFICIAL_MODE + ".jsonl"
            ),
        wildcard_constraints:
            evaluator="truvari|aardvark|vcfdist"
        threads:
            4
        resources:
            mem_mb=8192
        params:
            output_dir=lambda wildcards: (
                f"results/evaluation/{wildcards.tool}/{wildcards.evaluator}"
            ),
        conda:
            "../envs/core.yaml"
        shell:
            """
            {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
              --rule-name evaluate_{wildcards.evaluator} \
              --job-key {SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE} \
              --run-id {RUN_ID:q} \
              --module-or-tool-id pgbench-core \
              --snakefile-path {input.rule_source:q} \
              --rule-source-path {input.rule_source:q} \
              --script-or-wrapper-path {input.script:q} \
              --config-snapshot {input.config:q} \
              --score-profile {input.score_profile:q} \
              --run-context {input.context:q} \
              --pangenome-manifest {input.pangenome_manifest:q} \
              --reference {input.reference:q} \
              --truth-profile {config[truth][primary]:q} \
              --snakemake-version {SNAKEMAKE_VERSION:q} \
              --execution-profile local \
              --random-seed {config[execution][random_seed]} \
              --threads {threads} \
              --resource mem_mb={resources.mem_mb} \
              --wildcard tool={wildcards.tool:q} \
              --wildcard evaluator={wildcards.evaluator:q} \
              --input {input.query:q} \
              --input {input.truth:q} \
              --input {input.truth_index:q} \
              --input {input.regions:q} \
              --input {input.reference:q} \
              --input {input.config:q} \
              --input {input.script:q} \
              --output {output.ledger:q} \
              --output {output.complete:q} \
              --upstream-manifest {input.link_rule_manifest:q} \
              --manifest-output {output.rule_manifest:q} \
              -- \
              {PYTHON_EXECUTABLE:q} {input.script:q} \
                --evaluator {wildcards.evaluator:q} \
                --config {input.config:q} \
                --query {input.query:q} \
                --truth {input.truth:q} \
                --reference {input.reference:q} \
                --regions {input.regions:q} \
                --output-dir {params.output_dir:q} \
                --threads {threads} \
              > {log:q} 2>&1
            """
