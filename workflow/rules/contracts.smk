def contract_tool_manifest(wildcards):
    for settings in EXTERNAL_SETTINGS:
        if settings["tool_id"] == wildcards.tool:
            return settings["tool_manifest"]
    raise WorkflowError(f"no tool manifest for {wildcards.tool!r}")


def contract_resolved_inputs(wildcards):
    return TOOL_SETTINGS_BY_ID[wildcards.tool]["resolved_inputs"]


def contract_tool_rule_manifest(wildcards):
    return TOOL_SETTINGS_BY_ID[wildcards.tool]["rule_manifest"]


rule validate_evaluator_semantics:
    input:
        cases="tests/fixtures/evaluator_semantics/cases.yaml",
        truvari="tests/fixtures/evaluator_semantics/truvari.tsv",
        aardvark="tests/fixtures/evaluator_semantics/aardvark.tsv",
        vcfdist="tests/fixtures/evaluator_semantics/vcfdist.tsv",
        evaluator_profile=config["catalogs"]["evaluator_profile"],
        context=RESULTS_ROOT + "/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        rule_source="workflow/rules/contracts.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/validate_evaluator_semantics.py",
        environment=CORE_ENV_SPEC,
        provenance_library=PROVENANCE_LIBRARY,
    output:
        validation=SEMANTIC_VALIDATION_JSON,
        rule_manifest=SEMANTIC_VALIDATION_RULE_MANIFEST,
    log:
        f"{LOG_ROOT}/rules/validate_evaluator_semantics/contract.log",
    benchmark:
        f"{BENCHMARK_ROOT}/rules/validate_evaluator_semantics/contract.jsonl",
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name validate_evaluator_semantics --job-key contract \
          --run-id {RUN_ID:q} --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} --run-context {input.context:q} \
          --reference {config[reference][fasta]:q} --truth-profile {config[truth][primary]:q} \
          --snakemake-version {SNAKEMAKE_VERSION:q} --execution-profile local \
          --random-seed {config[execution][random_seed]} --threads 1 --resource mem_mb=1024 \
          --input {input.cases:q} --input {input.truvari:q} --input {input.aardvark:q} \
          --input {input.vcfdist:q} --input {input.evaluator_profile:q} --input {input.context:q} \
          --input {input.config:q} --input {input.score_profile:q} --input {input.rule_source:q} \
          --input {input.rule_executor:q} --input {input.script:q} --input {input.environment:q} \
          --input {input.provenance_library:q} --output {output.validation:q} \
          --upstream-manifest {VALIDATE_MANIFEST:q} --manifest-output {output.rule_manifest:q} -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} --cases {input.cases:q} \
            --truvari-ledger {input.truvari:q} --aardvark-ledger {input.aardvark:q} \
            --vcfdist-ledger {input.vcfdist:q} --output {output.validation:q} \
          > {log:q} 2>&1
        """


rule build_coverage_datasets:
    input:
        r1=config["sample"]["fastq_r1"],
        r2=config["sample"]["fastq_r2"],
        context=RESULTS_ROOT + "/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        rule_source="workflow/rules/contracts.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/build_coverage_fastqs.py",
        environment=CORE_ENV_SPEC,
        provenance_library=PROVENANCE_LIBRARY,
    output:
        subsets=directory(RESULTS_ROOT + "/coverage/subsets"),
        manifest=COVERAGE_MANIFEST,
        rule_manifest=COVERAGE_RULE_MANIFEST,
    log:
        f"{LOG_ROOT}/rules/build_coverage_datasets/coverage.log",
    benchmark:
        f"{BENCHMARK_ROOT}/rules/build_coverage_datasets/coverage.jsonl",
    params:
        coverage_args=lambda wildcards: cli_repeated("--coverage", [10, 20, 30]),
        seed_args=lambda wildcards: cli_repeated(
            "--seed", config["benchmark_contract"]["downsampling_seeds"]
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name build_coverage_datasets --job-key coverage \
          --run-id {RUN_ID:q} --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} --run-context {input.context:q} \
          --reference {config[reference][fasta]:q} --truth-profile {config[truth][primary]:q} \
          --snakemake-version {SNAKEMAKE_VERSION:q} --execution-profile local \
          --random-seed {config[execution][random_seed]} --threads 1 --resource mem_mb=4096 \
          --input {input.r1:q} --input {input.r2:q} --input {input.context:q} \
          --input {input.config:q} --input {input.score_profile:q} --input {input.rule_source:q} \
          --input {input.rule_executor:q} --input {input.script:q} --input {input.environment:q} \
          --input {input.provenance_library:q} --output {output.subsets:q} \
          --output {output.manifest:q} --upstream-manifest {CONTEXT_MANIFEST:q} \
          --manifest-output {output.rule_manifest:q} -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} --fastq-r1 {input.r1:q} --fastq-r2 {input.r2:q} \
            --output-dir {output.subsets:q} --source-coverage {config[benchmark_contract][source_coverage_x]} \
            {params.coverage_args:q} {params.seed_args:q} --manifest {output.manifest:q} \
          > {log:q} 2>&1
        """


rule freeze_information_contract:
    input:
        tool_manifest=contract_tool_manifest,
        resolved=contract_resolved_inputs,
        tool_rule_manifest=contract_tool_rule_manifest,
        context=RESULTS_ROOT + "/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        allowed_information_policy=config["catalogs"]["allowed_information"],
        tuning_policy=config["catalogs"]["tuning_policy"],
        rule_source="workflow/rules/contracts.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/freeze_information_contract.py",
        environment=CORE_ENV_SPEC,
        provenance_library=PROVENANCE_LIBRARY,
    output:
        allowed=RESULTS_ROOT + "/provenance/{tool}/contracts/allowed_inputs.json",
        parameters=RESULTS_ROOT + "/provenance/{tool}/contracts/parameter_manifest.json",
        resources=RESULTS_ROOT + "/provenance/{tool}/contracts/external_resource_hashes.json",
        tuning=RESULTS_ROOT + "/provenance/{tool}/contracts/training_or_tuning_status.json",
        rule_manifest=RESULTS_ROOT + "/provenance/rules/freeze_information_contract/{tool}.json",
    log:
        f"{LOG_ROOT}/rules/freeze_information_contract/" + "{tool}.log",
    benchmark:
        f"{BENCHMARK_ROOT}/rules/freeze_information_contract/" + "{tool}.jsonl",
    params:
        output_dir=RESULTS_ROOT + "/provenance/{tool}/contracts",
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name freeze_information_contract --job-key {wildcards.tool} \
          --run-id {RUN_ID:q} --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} --run-context {input.context:q} \
          --reference {config[reference][fasta]:q} --truth-profile {config[truth][primary]:q} \
          --snakemake-version {SNAKEMAKE_VERSION:q} --execution-profile local \
          --random-seed {config[execution][random_seed]} --threads 1 --resource mem_mb=1024 \
          --input {input.tool_manifest:q} --input {input.resolved:q} --input {input.context:q} \
          --input {input.config:q} --input {input.score_profile:q} \
          --input {input.allowed_information_policy:q} --input {input.tuning_policy:q} \
          --input {input.rule_source:q} \
          --input {input.rule_executor:q} --input {input.script:q} --input {input.environment:q} \
          --input {input.provenance_library:q} --output {output.allowed:q} \
          --output {output.parameters:q} --output {output.resources:q} --output {output.tuning:q} \
          --upstream-manifest {input.tool_rule_manifest:q} --manifest-output {output.rule_manifest:q} -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} --tool-manifest {input.tool_manifest:q} \
            --resolved-inputs {input.resolved:q} \
            --allowed-information-policy {input.allowed_information_policy:q} \
            --tuning-policy {input.tuning_policy:q} --output-dir {params.output_dir:q} \
          > {log:q} 2>&1
        """
