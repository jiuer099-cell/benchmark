def synthetic_metric_fixture(wildcards):
    if not SYNTHETIC_MODE:
        return []
    directory = config.get("development", {}).get("evaluator_fixture_dir")
    if not directory:
        raise WorkflowError(
            "synthetic scoring requires development.evaluator_fixture_dir"
        )
    return f"{directory}/score_input.json"


def formal_vote_ledgers(wildcards):
    if SYNTHETIC_MODE:
        return []
    return [
        f"results/evaluation/{wildcards.tool}/{evaluator}/votes.tsv"
        for evaluator in FORMAL_EVALUATORS
    ]


def formal_evaluator_manifests(wildcards):
    if SYNTHETIC_MODE:
        return []
    return [
        semantic_rule_manifest(
            f"evaluate_{evaluator}",
            f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}",
        )
        for evaluator in FORMAL_EVALUATORS
    ]


def metrics_materializer(wildcards):
    if SYNTHETIC_MODE:
        return "workflow/scripts/materialize_synthetic_metrics.py"
    return "workflow/scripts/materialize_formal_consensus_metrics.py"


def dynamic_metric_inputs(wildcards):
    return [
        *(
            [synthetic_metric_fixture(wildcards)]
            if SYNTHETIC_MODE
            else formal_vote_ledgers(wildcards)
        ),
        *formal_evaluator_manifests(wildcards),
    ]


def dynamic_metric_upstreams(wildcards):
    return [
        semantic_rule_manifest(
            "link_pangenome_alleles",
            f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}",
        ),
        *formal_evaluator_manifests(wildcards),
    ]


def materializer_arguments(wildcards):
    output = f"results/summary/{wildcards.tool}/metrics.json"
    if SYNTHETIC_MODE:
        return [
            "--fixture",
            synthetic_metric_fixture(wildcards),
            "--metric-dictionary",
            config["catalogs"]["metric_dictionary"],
            "--metrics-schema",
            "workflow/schemas/metrics.schema.yaml",
            "--score-profile",
            config["catalogs"]["score_weights"],
            "--evaluator-manifest",
            (
                "results/provenance/rules/"
                f"tool__{wildcards.tool}__execute/"
                f"{SAMPLE_ID}.{OFFICIAL_MODE}.json"
            ),
            "--pangenome-manifest",
            (
                "results/provenance/rules/link_pangenome_alleles/"
                f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}.json"
            ),
            "--resource-manifest",
            (
                "results/provenance/rules/"
                f"tool__{wildcards.tool}__execute/"
                f"{SAMPLE_ID}.{OFFICIAL_MODE}.json"
            ),
            "--output",
            output,
        ]
    evaluator_root = f"results/evaluation/{wildcards.tool}"
    return [
        "--truvari-ledger",
        f"{evaluator_root}/truvari/votes.tsv",
        "--aardvark-ledger",
        f"{evaluator_root}/aardvark/votes.tsv",
        "--vcfdist-ledger",
        f"{evaluator_root}/vcfdist/votes.tsv",
        *cli_repeated(
            "--evaluator-manifest",
            formal_evaluator_manifests(wildcards),
        ),
        "--score-profile",
        config["catalogs"]["score_weights"],
        "--run-id",
        RUN_ID,
        "--sample-id",
        SAMPLE_ID,
        "--tool-id",
        wildcards.tool,
        "--official-score-mode",
        OFFICIAL_MODE,
        "--primary-truth-profile",
        config["truth"]["primary"],
        "--output",
        output,
    ]


def fused_metrics_rule_manifest(wildcards):
    return semantic_rule_manifest(
        "fuse_evaluator_metrics",
        f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}",
    )


rule fuse_evaluator_metrics:
    input:
        fixture=synthetic_metric_fixture,
        evaluator_ledgers=formal_vote_ledgers,
        evaluator_manifests=formal_evaluator_manifests,
        linked=(
            f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            "{tool}/canonical/linked.vcf"
        ),
        links=(
            f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            "{tool}/canonical/allele_links.tsv"
        ),
        link_rule_manifest=(
            "results/provenance/rules/link_pangenome_alleles/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
        tool_rule_manifest=(
            "results/provenance/rules/tool__{tool}__execute/"
            f"{SAMPLE_ID}.{OFFICIAL_MODE}.json"
        ),
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        metric_dictionary=config["catalogs"]["metric_dictionary"],
        metrics_schema="workflow/schemas/metrics.schema.yaml",
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_executor=RULE_EXECUTOR,
        rule_source="workflow/rules/scoring.smk",
        core_env="workflow/envs/core.yaml",
        script=metrics_materializer,
        provenance_library=PROVENANCE_LIBRARY,
        metrics_library=METRICS_LIBRARY,
        scoring_library=SCORING_LIBRARY,
    output:
        metrics="results/summary/{tool}/metrics.json",
        rule_manifest=(
            "results/provenance/rules/fuse_evaluator_metrics/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
    log:
        (
            f"logs/rules/fuse_evaluator_metrics/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".log"
        ),
    benchmark:
        (
            f"benchmarks/rules/fuse_evaluator_metrics/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".jsonl"
        ),
    params:
        dynamic_input_args=lambda wildcards: cli_repeated(
            "--input", dynamic_metric_inputs(wildcards)
        ),
        dynamic_upstream_args=lambda wildcards: cli_repeated(
            "--upstream-manifest", dynamic_metric_upstreams(wildcards)
        ),
        materializer_args=materializer_arguments,
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name fuse_evaluator_metrics \
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
          --threads 1 \
          --resource mem_mb=512 \
          --wildcard tool={wildcards.tool:q} \
          --param evaluation_mode={EVALUATION_MODE:q} \
          {params.dynamic_input_args:q} \
          --input {input.linked:q} \
          --input {input.links:q} \
          --input {input.link_rule_manifest:q} \
          --input {input.tool_rule_manifest:q} \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.score_profile:q} \
          --input {input.metric_dictionary:q} \
          --input {input.metrics_schema:q} \
          --input {input.pangenome_manifest:q} \
          --input {input.reference:q} \
          --input {input.rule_executor:q} \
          --input {input.rule_source:q} \
          --input {input.core_env:q} \
          --input {input.script:q} \
          --input {input.provenance_library:q} \
          --input {input.metrics_library:q} \
          --input {input.scoring_library:q} \
          --output {output.metrics:q} \
          {params.dynamic_upstream_args:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            {params.materializer_args:q} \
          > {log:q} 2>&1
        """


rule compute_pgbench_score:
    input:
        metrics="results/summary/{tool}/metrics.json",
        metrics_rule_manifest=fused_metrics_rule_manifest,
        audit="results/provenance/{tool}/pre-score-audit.json",
        audit_rule_manifest=(
            "results/provenance/rules/audit_score_inputs/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        score_profile=config["catalogs"]["score_weights"],
        metric_dictionary=config["catalogs"]["metric_dictionary"],
        metrics_schema="workflow/schemas/metrics.schema.yaml",
        rule_executor=RULE_EXECUTOR,
        rule_source="workflow/rules/scoring.smk",
        core_env="workflow/envs/core.yaml",
        script="workflow/scripts/score_tools.py",
        provenance_library=PROVENANCE_LIBRARY,
        metrics_library=METRICS_LIBRARY,
        scoring_library=SCORING_LIBRARY,
    output:
        score="results/summary/{tool}/score.json",
        rule_manifest=(
            "results/provenance/rules/compute_pgbench_score/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
    log:
        (
            f"logs/rules/compute_pgbench_score/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".log"
        ),
    benchmark:
        (
            f"benchmarks/rules/compute_pgbench_score/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".jsonl"
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name compute_pgbench_score \
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
          --threads 1 \
          --resource mem_mb=512 \
          --wildcard tool={wildcards.tool:q} \
          --param score_profile={input.score_profile:q} \
          --param evaluation_mode={EVALUATION_MODE:q} \
          --input {input.metrics:q} \
          --input {input.audit:q} \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.pangenome_manifest:q} \
          --input {input.reference:q} \
          --input {input.score_profile:q} \
          --input {input.metric_dictionary:q} \
          --input {input.metrics_schema:q} \
          --input {input.rule_executor:q} \
          --input {input.rule_source:q} \
          --input {input.core_env:q} \
          --input {input.script:q} \
          --input {input.provenance_library:q} \
          --input {input.metrics_library:q} \
          --input {input.scoring_library:q} \
          --output {output.score:q} \
          --upstream-manifest {input.metrics_rule_manifest:q} \
          --upstream-manifest {input.audit_rule_manifest:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            --metrics {input.metrics:q} \
            --provenance-audit {input.audit:q} \
            --run-context {input.context:q} \
            --score-profile {input.score_profile:q} \
            --metric-dictionary {input.metric_dictionary:q} \
            --metrics-schema {input.metrics_schema:q} \
            --evaluation-mode {EVALUATION_MODE:q} \
            --expected-run-id {RUN_ID:q} \
            --expected-sample-id {SAMPLE_ID:q} \
            --expected-tool-id {wildcards.tool:q} \
            --expected-official-score-mode {OFFICIAL_MODE:q} \
            --expected-primary-truth-profile {config[truth][primary]:q} \
            --output {output.score:q} \
          > {log:q} 2>&1
        """
