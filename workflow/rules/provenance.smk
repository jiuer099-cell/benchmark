def pre_score_manifest_paths(wildcards):
    tool_job = f"{SAMPLE_ID}.{OFFICIAL_MODE}"
    core_job = f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}"
    return [
        VALIDATE_MANIFEST,
        CONTEXT_MANIFEST,
        PANGENOME_RULE_MANIFEST,
        CHALLENGE_RULE_MANIFEST,
        semantic_rule_manifest(
            f"tool__{wildcards.tool}__execute",
            tool_job,
        ),
        semantic_rule_manifest("canonicalize_vcf", core_job),
        semantic_rule_manifest("link_pangenome_alleles", core_job),
        semantic_rule_manifest("fuse_evaluator_metrics", core_job),
    ]


def expected_pre_score_jobs(wildcards):
    tool_job = f"{SAMPLE_ID}.{OFFICIAL_MODE}"
    core_job = f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}"
    return [
        "validate_config=config",
        "snapshot_run_context=context",
        f"build_pangenome_manifest={PANGENOME_ID}",
        f"build_blinded_challenge_panel={SAMPLE_ID}",
        f"tool__{wildcards.tool}__execute={tool_job}",
        f"canonicalize_vcf={core_job}",
        f"link_pangenome_alleles={core_job}",
        f"fuse_evaluator_metrics={core_job}",
    ]


def pre_score_artifact_paths(wildcards):
    tool_root = f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/{wildcards.tool}"
    return [
        "results/provenance/config.validated.json",
        "results/provenance/run-context.json",
        f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        (
            f"results/pangenome/{PANGENOME_ID}/challenge/"
            f"{SAMPLE_ID}.blinded.vcf"
        ),
        external_raw_vcf(wildcards),
        f"{tool_root}/canonical/calls.vcf",
        f"{tool_root}/canonical/linked.vcf",
        f"results/summary/{wildcards.tool}/metrics.json",
    ]


rule build_rule_lineage:
    input:
        manifests=pre_score_manifest_paths,
        artifacts=pre_score_artifact_paths,
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_executor=RULE_EXECUTOR,
        rule_source="workflow/rules/provenance.smk",
        core_env="workflow/envs/core.yaml",
        script="workflow/scripts/build_rule_lineage.py",
        provenance_library=PROVENANCE_LIBRARY,
    output:
        json="results/provenance/{tool}/pre-score-lineage.json",
        tsv="results/provenance/{tool}/pre-score-lineage.tsv",
        rule_manifest=(
            "results/provenance/rules/build_rule_lineage/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
    log:
        (
            f"logs/rules/build_rule_lineage/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".log"
        ),
    benchmark:
        (
            f"benchmarks/rules/build_rule_lineage/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".jsonl"
        ),
    params:
        manifest_args=lambda wildcards: cli_repeated(
            "--manifest", pre_score_manifest_paths(wildcards)
        ),
        provenance_input_args=lambda wildcards: cli_repeated(
            "--input",
            [
                *pre_score_manifest_paths(wildcards),
                *pre_score_artifact_paths(wildcards),
                "results/provenance/run-context.json",
                CONFIG_PATH,
                config["catalogs"]["score_weights"],
                f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
                config["reference"]["fasta"],
                RULE_EXECUTOR,
                "workflow/rules/provenance.smk",
                "workflow/envs/core.yaml",
                "workflow/scripts/build_rule_lineage.py",
                PROVENANCE_LIBRARY,
            ],
        ),
        upstream_args=lambda wildcards: cli_repeated(
            "--upstream-manifest", pre_score_manifest_paths(wildcards)
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name build_rule_lineage \
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
          {params.provenance_input_args:q} \
          --output {output.json:q} \
          --output {output.tsv:q} \
          {params.upstream_args:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            {params.manifest_args:q} \
            --output-json {output.json:q} \
            --output-tsv {output.tsv:q} \
          > {log:q} 2>&1
        """


rule audit_score_inputs:
    input:
        manifests=pre_score_manifest_paths,
        artifacts=pre_score_artifact_paths,
        lineage="results/provenance/{tool}/pre-score-lineage.json",
        lineage_rule_manifest=(
            "results/provenance/rules/build_rule_lineage/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_executor=RULE_EXECUTOR,
        rule_source="workflow/rules/provenance.smk",
        core_env="workflow/envs/core.yaml",
        script="workflow/scripts/audit_provenance.py",
        provenance_library=PROVENANCE_LIBRARY,
    output:
        audit="results/provenance/{tool}/pre-score-audit.json",
        rule_manifest=(
            "results/provenance/rules/audit_score_inputs/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
    log:
        (
            f"logs/rules/audit_score_inputs/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".log"
        ),
    benchmark:
        (
            f"benchmarks/rules/audit_score_inputs/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".jsonl"
        ),
    params:
        manifest_args=lambda wildcards: cli_repeated(
            "--manifest", pre_score_manifest_paths(wildcards)
        ),
        expected_args=lambda wildcards: cli_repeated(
            "--expected-job", expected_pre_score_jobs(wildcards)
        ),
        provenance_input_args=lambda wildcards: cli_repeated(
            "--input",
            [
                *pre_score_manifest_paths(wildcards),
                *pre_score_artifact_paths(wildcards),
                f"results/provenance/{wildcards.tool}/pre-score-lineage.json",
                "results/provenance/run-context.json",
                CONFIG_PATH,
                config["catalogs"]["score_weights"],
                f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
                config["reference"]["fasta"],
                RULE_EXECUTOR,
                "workflow/rules/provenance.smk",
                "workflow/envs/core.yaml",
                "workflow/scripts/audit_provenance.py",
                PROVENANCE_LIBRARY,
            ],
        ),
        upstream_args=lambda wildcards: cli_repeated(
            "--upstream-manifest",
            [
                *pre_score_manifest_paths(wildcards),
                semantic_rule_manifest(
                    "build_rule_lineage",
                    f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}",
                ),
            ],
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name audit_score_inputs \
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
          {params.provenance_input_args:q} \
          --output {output.audit:q} \
          {params.upstream_args:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            {params.manifest_args:q} \
            {params.expected_args:q} \
            --workspace-root . \
            --verify-paths \
            --output {output.audit:q} \
          > {log:q} 2>&1
        """


def final_score_supporting_manifest_paths(wildcards):
    return [
        *pre_score_manifest_paths(wildcards),
        semantic_rule_manifest(
            "build_rule_lineage",
            f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}",
        ),
    ]


rule finalize_score_provenance:
    input:
        score="results/summary/{tool}/score.json",
        metrics="results/summary/{tool}/metrics.json",
        score_rule_manifest=(
            "results/provenance/rules/compute_pgbench_score/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
        audit="results/provenance/{tool}/pre-score-audit.json",
        audit_rule_manifest=(
            "results/provenance/rules/audit_score_inputs/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
        supporting_manifests=final_score_supporting_manifest_paths,
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_executor=RULE_EXECUTOR,
        rule_source="workflow/rules/provenance.smk",
        core_env="workflow/envs/core.yaml",
        script="workflow/scripts/finalize_score_provenance.py",
        provenance_library=PROVENANCE_LIBRARY,
    output:
        package="results/summary/{tool}/score-package.json",
        lineage_json="results/provenance/{tool}/rule-lineage.json",
        lineage_tsv="results/provenance/{tool}/rule-lineage.tsv",
        audit="results/provenance/{tool}/provenance-audit.json",
        rule_manifest=(
            "results/provenance/rules/finalize_score_provenance/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
    log:
        (
            f"logs/rules/finalize_score_provenance/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".log"
        ),
    benchmark:
        (
            f"benchmarks/rules/finalize_score_provenance/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".jsonl"
        ),
    params:
        manifest_args=lambda wildcards: cli_repeated(
            "--manifest", final_score_supporting_manifest_paths(wildcards)
        ),
        provenance_input_args=lambda wildcards: cli_repeated(
            "--input",
            [
                f"results/summary/{wildcards.tool}/score.json",
                f"results/summary/{wildcards.tool}/metrics.json",
                semantic_rule_manifest(
                    "compute_pgbench_score",
                    f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}",
                ),
                f"results/provenance/{wildcards.tool}/pre-score-audit.json",
                semantic_rule_manifest(
                    "audit_score_inputs",
                    f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}",
                ),
                *final_score_supporting_manifest_paths(wildcards),
                "results/provenance/run-context.json",
                CONFIG_PATH,
                config["catalogs"]["score_weights"],
                f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
                config["reference"]["fasta"],
                RULE_EXECUTOR,
                "workflow/rules/provenance.smk",
                "workflow/envs/core.yaml",
                "workflow/scripts/finalize_score_provenance.py",
                PROVENANCE_LIBRARY,
            ],
        ),
        upstream_args=lambda wildcards: cli_repeated(
            "--upstream-manifest",
            [
                semantic_rule_manifest(
                    "audit_score_inputs",
                    f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}",
                ),
                semantic_rule_manifest(
                    "compute_pgbench_score",
                    f"{SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE}",
                ),
            ],
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name finalize_score_provenance \
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
          {params.provenance_input_args:q} \
          --output {output.package:q} \
          --output {output.lineage_json:q} \
          --output {output.lineage_tsv:q} \
          --output {output.audit:q} \
          {params.upstream_args:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            --score-json {input.score:q} \
            --score-rule-manifest {input.score_rule_manifest:q} \
            --pre-score-audit {input.audit:q} \
            --audit-rule-manifest {input.audit_rule_manifest:q} \
            {params.manifest_args:q} \
            --workspace-root . \
            --output-package {output.package:q} \
            --output-lineage-json {output.lineage_json:q} \
            --output-lineage-tsv {output.lineage_tsv:q} \
            --output-audit {output.audit:q} \
          > {log:q} 2>&1
        """
