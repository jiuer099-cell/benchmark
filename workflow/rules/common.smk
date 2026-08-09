rule validate_config:
    input:
        config=CONFIG_PATH,
        config_schema="config/config.schema.yaml",
        tool_schema="workflow/schemas/tool.schema.yaml",
        rule_source="workflow/rules/common.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/validate_config.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    params:
        score_profile=config["catalogs"]["score_weights"],
    output:
        data=RESULTS_ROOT + "/provenance/config.validated.json",
        rule_manifest=VALIDATE_MANIFEST,
    log:
        LOG_ROOT + "/rules/validate_config/config.log",
    benchmark:
        BENCHMARK_ROOT + "/rules/validate_config/config.jsonl",
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name validate_config \
          --job-key config \
          --run-id {RUN_ID:q} \
          --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} \
          --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} \
          --config-snapshot {input.config:q} \
          --score-profile {params.score_profile:q} \
          --truth-profile {config[truth][primary]:q} \
          --snakemake-version {SNAKEMAKE_VERSION:q} \
          --execution-profile local \
          --random-seed {config[execution][random_seed]} \
          --threads 1 \
          --resource mem_mb=512 \
          --param config_schema={input.config_schema:q} \
          --input {input.config:q} \
          --input {input.config_schema:q} \
          --input {input.tool_schema:q} \
          --input {input.rule_source:q} \
          --input {input.rule_executor:q} \
          --input {input.script:q} \
          --input {input.environment:q} \
          --input {input.provenance_library:q} \
          --output {output.data:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            --config {input.config:q} \
            --config-schema {input.config_schema:q} \
            --tool-schema {input.tool_schema:q} \
            --repo-root . \
            --output {output.data:q} \
          > {log:q} 2>&1
        """


rule snapshot_run_context:
    input:
        validated=RESULTS_ROOT + "/provenance/config.validated.json",
        validated_manifest=VALIDATE_MANIFEST,
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        evaluator_profile=config["catalogs"]["evaluator_profile"],
        rule_source="workflow/rules/common.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/snapshot_run_context.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    output:
        data=RESULTS_ROOT + "/provenance/run-context.json",
        rule_manifest=CONTEXT_MANIFEST,
    log:
        LOG_ROOT + "/rules/snapshot_run_context/context.log",
    benchmark:
        BENCHMARK_ROOT + "/rules/snapshot_run_context/context.jsonl",
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name snapshot_run_context \
          --job-key context \
          --run-id {RUN_ID:q} \
          --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} \
          --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} \
          --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} \
          --run-context {output.data:q} \
          --truth-profile {config[truth][primary]:q} \
          --snakemake-version {SNAKEMAKE_VERSION:q} \
          --execution-profile local \
          --random-seed {config[execution][random_seed]} \
          --threads 1 \
          --resource mem_mb=512 \
          --param score_profile={input.score_profile:q} \
          --param evaluator_profile={input.evaluator_profile:q} \
          --input {input.validated:q} \
          --input {input.config:q} \
          --input {input.score_profile:q} \
          --input {input.evaluator_profile:q} \
          --input {input.rule_source:q} \
          --input {input.rule_executor:q} \
          --input {input.script:q} \
          --input {input.environment:q} \
          --input {input.provenance_library:q} \
          --output {output.data:q} \
          --upstream-manifest {input.validated_manifest:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            --repo-root . \
            --score-profile {input.score_profile:q} \
            --evaluator-profile {input.evaluator_profile:q} \
            --random-seed {config[execution][random_seed]} \
            --snakemake-version {SNAKEMAKE_VERSION:q} \
            --execution-profile local \
            --output {output.data:q} \
          > {log:q} 2>&1
        """
