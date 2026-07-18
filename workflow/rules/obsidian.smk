OBSIDIAN_JOB_KEY = config["project"]["id"]


rule sync_obsidian_design:
    input:
        context="results/provenance/run-context.json",
        context_rule_manifest=CONTEXT_MANIFEST,
        config=CONFIG_PATH,
        design=(
            "docs/superpowers/specs/"
            "2026-07-16-hg002-grch37-sv-benchmark-design.md"
        ),
        binding=config["obsidian"]["binding_file"],
        rule_registry="workflow/rule-registry.yaml",
        score_profile=config["catalogs"]["score_weights"],
        pangenome_schema=(
            "workflow/schemas/pangenome-manifest.schema.yaml"
        ),
        rule_source="workflow/rules/obsidian.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/sync_obsidian_design.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    output:
        check="results/provenance/obsidian-sync.json",
        rule_manifest=(
            "results/provenance/rules/sync_obsidian_design/"
            f"{OBSIDIAN_JOB_KEY}.json"
        ),
    log:
        f"logs/rules/sync_obsidian_design/{OBSIDIAN_JOB_KEY}.log",
    benchmark:
        (
            "benchmarks/rules/sync_obsidian_design/"
            f"{OBSIDIAN_JOB_KEY}.jsonl"
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        python {input.rule_executor:q} \
          --rule-name sync_obsidian_design \
          --job-key {OBSIDIAN_JOB_KEY:q} \
          --run-id {RUN_ID:q} \
          --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} \
          --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} \
          --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} \
          --run-context {input.context:q} \
          --truth-profile {config[truth][primary]:q} \
          --snakemake-version {SNAKEMAKE_VERSION:q} \
          --execution-profile local \
          --random-seed {config[execution][random_seed]} \
          --threads 1 \
          --resource mem_mb=512 \
          --param code_state=partial \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.design:q} \
          --input {input.binding:q} \
          --input {input.rule_registry:q} \
          --input {input.score_profile:q} \
          --input {input.pangenome_schema:q} \
          --input {input.rule_source:q} \
          --input {input.rule_executor:q} \
          --input {input.script:q} \
          --input {input.environment:q} \
          --input {input.provenance_library:q} \
          --output {output.check:q} \
          --upstream-manifest {input.context_rule_manifest:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          python {input.script:q} \
            --repo-root . \
            --design {input.design:q} \
            --config {input.config:q} \
            --rule-registry {input.rule_registry:q} \
            --binding {input.binding:q} \
            --pangenome-schema {input.pangenome_schema:q} \
            --code-state partial \
            --output-check {output.check:q} \
          > {log:q} 2>&1
        """


rule check_obsidian_sync:
    input:
        context="results/provenance/run-context.json",
        context_rule_manifest=CONTEXT_MANIFEST,
        config=CONFIG_PATH,
        design=(
            "docs/superpowers/specs/"
            "2026-07-16-hg002-grch37-sv-benchmark-design.md"
        ),
        binding=config["obsidian"]["binding_file"],
        rule_registry="workflow/rule-registry.yaml",
        score_profile=config["catalogs"]["score_weights"],
        pangenome_schema=(
            "workflow/schemas/pangenome-manifest.schema.yaml"
        ),
        rule_source="workflow/rules/obsidian.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/check_obsidian_sync.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    output:
        check="results/provenance/obsidian-sync-check.json",
        rule_manifest=(
            "results/provenance/rules/check_obsidian_sync/"
            f"{OBSIDIAN_JOB_KEY}.json"
        ),
    log:
        f"logs/rules/check_obsidian_sync/{OBSIDIAN_JOB_KEY}.log",
    benchmark:
        (
            "benchmarks/rules/check_obsidian_sync/"
            f"{OBSIDIAN_JOB_KEY}.jsonl"
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        python {input.rule_executor:q} \
          --rule-name check_obsidian_sync \
          --job-key {OBSIDIAN_JOB_KEY:q} \
          --run-id {RUN_ID:q} \
          --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} \
          --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} \
          --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} \
          --run-context {input.context:q} \
          --truth-profile {config[truth][primary]:q} \
          --snakemake-version {SNAKEMAKE_VERSION:q} \
          --execution-profile local \
          --random-seed {config[execution][random_seed]} \
          --threads 1 \
          --resource mem_mb=512 \
          --param code_state=partial \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.design:q} \
          --input {input.binding:q} \
          --input {input.rule_registry:q} \
          --input {input.score_profile:q} \
          --input {input.pangenome_schema:q} \
          --input {input.rule_source:q} \
          --input {input.rule_executor:q} \
          --input {input.script:q} \
          --input {input.environment:q} \
          --input {input.provenance_library:q} \
          --output {output.check:q} \
          --upstream-manifest {input.context_rule_manifest:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          python {input.script:q} \
            --repo-root . \
            --design {input.design:q} \
            --config {input.config:q} \
            --rule-registry {input.rule_registry:q} \
            --binding {input.binding:q} \
            --pangenome-schema {input.pangenome_schema:q} \
            --code-state partial \
            --output {output.check:q} \
          > {log:q} 2>&1
        """
