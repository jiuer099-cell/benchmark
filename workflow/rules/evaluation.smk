FORMAL_EVALUATORS = ("truvari", "aardvark", "vcfdist")


def evaluation_tool_manifest(wildcards):
    for settings in EXTERNAL_SETTINGS:
        if settings["tool_id"] == wildcards.tool:
            return settings["tool_manifest"]
    raise WorkflowError(f"no tool manifest for {wildcards.tool!r}")


if not SYNTHETIC_MODE:
    rule prepare_primary_truth:
        input:
            source=config["evaluation"]["truth_vcf"],
            source_index=config["evaluation"]["truth_vcf"] + ".tbi",
            regions=config["evaluation"]["benchmark_bed"],
            catalog=config["catalogs"]["truthsets"],
            evaluator_profile=config["catalogs"]["evaluator_profile"],
            context=RESULTS_ROOT + "/provenance/run-context.json",
            config=CONFIG_PATH,
            score_profile=config["catalogs"]["score_weights"],
            reference=config["reference"]["fasta"],
            rule_source="workflow/rules/evaluation.smk",
            rule_executor=RULE_EXECUTOR,
            script="workflow/scripts/prepare_sv_truth.py",
            environment=CORE_ENV_SPEC,
            provenance_library=PROVENANCE_LIBRARY,
        output:
            truth=PRIMARY_TRUTH_VCF,
            index=PRIMARY_TRUTH_INDEX,
            audit=PRIMARY_TRUTH_AUDIT,
            rule_manifest=PRIMARY_TRUTH_RULE_MANIFEST,
        log:
            (
                f"{LOG_ROOT}/rules/prepare_primary_truth/"
                f"{config['truth']['primary']}.log"
            ),
        benchmark:
            (
                f"{BENCHMARK_ROOT}/rules/prepare_primary_truth/"
                f"{config['truth']['primary']}.jsonl"
            ),
        conda:
            "../envs/core.yaml"
        shell:
            """
            {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
              --rule-name prepare_primary_truth \
              --job-key {config[truth][primary]:q} \
              --run-id {RUN_ID:q} \
              --module-or-tool-id pgbench-core \
              --snakefile-path {input.rule_source:q} \
              --rule-source-path {input.rule_source:q} \
              --script-or-wrapper-path {input.script:q} \
              --config-snapshot {input.config:q} \
              --score-profile {input.score_profile:q} \
              --run-context {input.context:q} \
              --reference {input.reference:q} \
              --truth-profile {config[truth][primary]:q} \
              --snakemake-version {SNAKEMAKE_VERSION:q} \
              --execution-profile local \
              --random-seed {config[execution][random_seed]} \
              --threads 1 \
              --resource mem_mb=4096 \
              --param minimum_sv_size=50 \
              --input {input.source:q} \
              --input {input.source_index:q} \
              --input {input.regions:q} \
              --input {input.catalog:q} \
              --input {input.evaluator_profile:q} \
              --input {input.context:q} \
              --input {input.config:q} \
              --input {input.score_profile:q} \
              --input {input.reference:q} \
              --input {input.rule_source:q} \
              --input {input.rule_executor:q} \
              --input {input.script:q} \
              --input {input.environment:q} \
              --input {input.provenance_library:q} \
              --output {output.truth:q} \
              --output {output.index:q} \
              --output {output.audit:q} \
              --upstream-manifest {CONTEXT_MANIFEST:q} \
              --manifest-output {output.rule_manifest:q} \
              -- \
              {PYTHON_EXECUTABLE:q} {input.script:q} \
                --source-vcf {input.source:q} \
                --source-index {input.source_index:q} \
                --benchmark-bed {input.regions:q} \
                --truth-catalog {input.catalog:q} \
                --evaluator-profile {input.evaluator_profile:q} \
                --profile-id {config[truth][primary]:q} \
                --output-vcf {output.truth:q} \
                --audit-json {output.audit:q} \
              > {log:q} 2>&1
            """

    if NOVEL_DISCOVERY_TOOL_IDS:
        rule prepare_novel_truth:
            input:
                truth=PRIMARY_TRUTH_VCF,
                truth_index=PRIMARY_TRUTH_INDEX,
                truth_rule_manifest=PRIMARY_TRUTH_RULE_MANIFEST,
                graph_gfa=config["pangenome"]["graph_assets"]["gfa"],
                variation_calls=(
                    config["pangenome"]["graph_assets"]["variation_calls"]
                ),
                graph_asset_lock=GRAPH_ASSET_LOCK,
                graph_rule_manifest=GRAPH_ASSET_RULE_MANIFEST,
                pangenome_manifest=(
                    f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml"
                ),
                pangenome_rule_manifest=PANGENOME_RULE_MANIFEST,
                novel_truth_profile=config["catalogs"]["novel_truth_profile"],
                evaluator_profile=config["catalogs"]["evaluator_profile"],
                reference=config["reference"]["fasta"],
                reference_index=config["reference"]["fai"],
                context=RESULTS_ROOT + "/provenance/run-context.json",
                config=CONFIG_PATH,
                score_profile=config["catalogs"]["score_weights"],
                rule_source="workflow/rules/evaluation.smk",
                rule_executor=RULE_EXECUTOR,
                script="workflow/scripts/prepare_novel_truth.py",
                environment=CORE_ENV_SPEC,
                provenance_library=PROVENANCE_LIBRARY,
            output:
                truth=NOVEL_TRUTH_VCF,
                index=NOVEL_TRUTH_INDEX,
                audit=NOVEL_TRUTH_AUDIT,
                exclusion_ledger=NOVEL_TRUTH_EXCLUSION_LEDGER,
                rule_manifest=NOVEL_TRUTH_RULE_MANIFEST,
            log:
                (
                    f"{LOG_ROOT}/rules/prepare_novel_truth/"
                    f"{config['truth']['primary']}.{PANGENOME_ID}.log"
                ),
            benchmark:
                (
                    f"{BENCHMARK_ROOT}/rules/prepare_novel_truth/"
                    f"{config['truth']['primary']}.{PANGENOME_ID}.jsonl"
                ),
            threads:
                1
            resources:
                mem_mb=16384
            conda:
                "../envs/core.yaml"
            shell:
                """
                {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
                  --rule-name prepare_novel_truth \
                  --job-key {config[truth][primary]}.{PANGENOME_ID} \
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
                  --truth-profile {NOVEL_TRUTH_PROFILE_ID:q} \
                  --snakemake-version {SNAKEMAKE_VERSION:q} \
                  --execution-profile local \
                  --random-seed {config[execution][random_seed]} \
                  --threads {threads} \
                  --resource mem_mb={resources.mem_mb} \
                  --param exclusion_policy=any_compatible_graph_allele \
                  --input {input.truth:q} \
                  --input {input.truth_index:q} \
                  --input {input.graph_gfa:q} \
                  --input {input.variation_calls:q} \
                  --input {input.graph_asset_lock:q} \
                  --input {input.pangenome_manifest:q} \
                  --input {input.novel_truth_profile:q} \
                  --input {input.evaluator_profile:q} \
                  --input {input.reference:q} \
                  --input {input.reference_index:q} \
                  --input {input.context:q} \
                  --input {input.config:q} \
                  --input {input.score_profile:q} \
                  --input {input.rule_source:q} \
                  --input {input.rule_executor:q} \
                  --input {input.script:q} \
                  --input {input.environment:q} \
                  --input {input.provenance_library:q} \
                  --output {output.truth:q} \
                  --output {output.index:q} \
                  --output {output.audit:q} \
                  --output {output.exclusion_ledger:q} \
                  --upstream-manifest {input.truth_rule_manifest:q} \
                  --upstream-manifest {input.graph_rule_manifest:q} \
                  --upstream-manifest {input.pangenome_rule_manifest:q} \
                  --manifest-output {output.rule_manifest:q} \
                  -- \
                  {PYTHON_EXECUTABLE:q} {input.script:q} \
                    --truth-vcf {input.truth:q} \
                    --reference {input.reference:q} \
                    --reference-index {input.reference_index:q} \
                    --graph-gfa {input.graph_gfa:q} \
                    --variation-calls {input.variation_calls:q} \
                    --graph-asset-lock {input.graph_asset_lock:q} \
                    --novel-truth-profile {input.novel_truth_profile:q} \
                    --evaluator-profile {input.evaluator_profile:q} \
                    --output-vcf {output.truth:q} \
                    --exclusion-ledger {output.exclusion_ledger:q} \
                    --audit-json {output.audit:q} \
                  > {log:q} 2>&1
                """

    rule run_formal_evaluator:
        input:
            query=(
                f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/"
                "{tool}/canonical/linked.vcf"
            ),
            link_rule_manifest=(
                RESULTS_ROOT + "/provenance/rules/link_pangenome_alleles/"
                f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
            ),
            truth=evaluation_truth_vcf,
            truth_index=evaluation_truth_index,
            truth_rule_manifest=evaluation_truth_rule_manifest,
            regions=config["evaluation"]["benchmark_bed"],
            reference=config["reference"]["fasta"],
            pangenome_manifest=(
                f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml"
            ),
            context=RESULTS_ROOT + "/provenance/run-context.json",
            config=CONFIG_PATH,
            score_profile=config["catalogs"]["score_weights"],
            evaluator_profile=config["catalogs"]["evaluator_profile"],
            tool_manifest=evaluation_tool_manifest,
            rule_source="workflow/rules/evaluation.smk",
            rule_executor=RULE_EXECUTOR,
            script="workflow/scripts/run_formal_evaluator.py",
            environment=CORE_ENV_SPEC,
            provenance_library=PROVENANCE_LIBRARY,
        output:
            ledger=RESULTS_ROOT + "/evaluation/{tool}/{evaluator}/votes.tsv",
            complete=RESULTS_ROOT + "/evaluation/{tool}/{evaluator}/.complete.json",
            rule_manifest=(
                RESULTS_ROOT + "/provenance/rules/evaluate_{evaluator}/"
                f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
            ),
        log:
            (
                f"{LOG_ROOT}/rules/evaluate_" + "{evaluator}/"
                f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".log"
            ),
        benchmark:
            (
                f"{BENCHMARK_ROOT}/rules/evaluate_" + "{evaluator}/"
                f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".jsonl"
            ),
        wildcard_constraints:
            evaluator="truvari|aardvark|vcfdist"
        threads:
            4
        resources:
            mem_mb=8192
        params:
            output_dir=lambda wildcards: (
                f"{RESULTS_ROOT}/evaluation/{wildcards.tool}/{wildcards.evaluator}"
            ),
            truth_profile=evaluation_truth_profile,
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
              --truth-profile {params.truth_profile:q} \
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
              --input {input.evaluator_profile:q} \
              --input {input.tool_manifest:q} \
              --input {input.script:q} \
              --output {output.ledger:q} \
              --output {output.complete:q} \
              --upstream-manifest {input.link_rule_manifest:q} \
              --upstream-manifest {input.truth_rule_manifest:q} \
              --manifest-output {output.rule_manifest:q} \
              -- \
              {PYTHON_EXECUTABLE:q} {input.script:q} \
                --evaluator {wildcards.evaluator:q} \
                --config {input.config:q} \
                --evaluator-profile {input.evaluator_profile:q} \
                --tool-manifest {input.tool_manifest:q} \
                --query {input.query:q} \
                --truth {input.truth:q} \
                --reference {input.reference:q} \
                --regions {input.regions:q} \
                --output-dir {params.output_dir:q} \
                --threads {threads} \
              > {log:q} 2>&1
            """
