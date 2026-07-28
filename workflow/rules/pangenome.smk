def population_vcf_input(wildcards):
    path = config.get("development", {}).get("population_vcf")
    if path:
        return path
    path = config.get("pangenome", {}).get("population_vcf")
    if path:
        return path
    raise WorkflowError("pangenome.population_vcf is required in formal mode")


def configured_graph_asset(name):
    def resolve(wildcards):
        value = config["pangenome"]["graph_assets"].get(name)
        if not value:
            raise WorkflowError(
                f"pangenome.graph_assets.{name} is required when "
                "build_graph_assets is true"
            )
        return value

    return resolve


rule lock_graph_assets:
    input:
        validated="results/provenance/config.validated.json",
        validated_manifest=VALIDATE_MANIFEST,
        context="results/provenance/run-context.json",
        context_manifest=CONTEXT_MANIFEST,
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        reference=config["reference"]["fasta"],
        source_manifest=configured_graph_asset("manifest"),
        gbz=configured_graph_asset("gbz"),
        xg=configured_graph_asset("xg"),
        min_index=configured_graph_asset("min"),
        dist=configured_graph_asset("dist"),
        sample_list=configured_graph_asset("sample_list"),
        rule_source="workflow/rules/pangenome.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/validate_graph_assets.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    output:
        lock=GRAPH_ASSET_LOCK,
        rule_manifest=GRAPH_ASSET_RULE_MANIFEST,
    log:
        f"logs/rules/lock_graph_assets/{PANGENOME_ID}.log",
    benchmark:
        f"benchmarks/rules/lock_graph_assets/{PANGENOME_ID}.jsonl",
    conda:
        "../envs/core.yaml"
    params:
        excluded_sample_csv=",".join(config["pangenome"]["excluded_samples"]),
        exclude_sample_arguments=cli_repeated(
            "--exclude-sample",
            config["pangenome"]["excluded_samples"],
        ),
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name lock_graph_assets \
          --job-key {PANGENOME_ID:q} \
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
          --resource mem_mb=1024 \
          --param reference_path={config[pangenome][graph_assets][reference_path]:q} \
          --param excluded_samples={params.excluded_sample_csv:q} \
          --input {input.validated:q} \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.score_profile:q} \
          --input {input.reference:q} \
          --input {input.source_manifest:q} \
          --input {input.gbz:q} \
          --input {input.xg:q} \
          --input {input.min_index:q} \
          --input {input.dist:q} \
          --input {input.sample_list:q} \
          --input {input.rule_source:q} \
          --input {input.rule_executor:q} \
          --input {input.script:q} \
          --input {input.environment:q} \
          --input {input.provenance_library:q} \
          --output {output.lock:q} \
          --upstream-manifest {input.validated_manifest:q} \
          --upstream-manifest {input.context_manifest:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            --source-manifest {input.source_manifest:q} \
            --gbz {input.gbz:q} \
            --xg {input.xg:q} \
            --min {input.min_index:q} \
            --dist {input.dist:q} \
            --sample-list {input.sample_list:q} \
            --reference-path {config[pangenome][graph_assets][reference_path]:q} \
            {params.exclude_sample_arguments:q} \
            --output {output.lock:q} \
          > {log:q} 2>&1
        """


rule build_pangenome_manifest:
    input:
        validated="results/provenance/config.validated.json",
        validated_manifest=VALIDATE_MANIFEST,
        context="results/provenance/run-context.json",
        context_manifest=CONTEXT_MANIFEST,
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        reference=config["reference"]["fasta"],
        population=population_vcf_input,
        graph_lock=([GRAPH_ASSET_LOCK] if GRAPH_ASSETS_ENABLED else []),
        graph_rule_manifest=(
            [GRAPH_ASSET_RULE_MANIFEST] if GRAPH_ASSETS_ENABLED else []
        ),
        rule_source="workflow/rules/pangenome.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/build_pangenome_manifest.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    output:
        panel=f"results/pangenome/{PANGENOME_ID}/panel.vcf",
        ledger=f"results/pangenome/{PANGENOME_ID}/allele-ledger.tsv",
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        rule_manifest=PANGENOME_RULE_MANIFEST,
    log:
        f"logs/rules/build_pangenome_manifest/{PANGENOME_ID}.log",
    benchmark:
        f"benchmarks/rules/build_pangenome_manifest/{PANGENOME_ID}.jsonl",
    conda:
        "../envs/core.yaml"
    params:
        excluded_sample_csv=",".join(config["pangenome"]["excluded_samples"]),
        graph_input_arguments=(
            ["--input", GRAPH_ASSET_LOCK] if GRAPH_ASSETS_ENABLED else []
        ),
        graph_upstream_arguments=(
            ["--upstream-manifest", GRAPH_ASSET_RULE_MANIFEST]
            if GRAPH_ASSETS_ENABLED
            else []
        ),
        graph_script_arguments=(
            ["--graph-assets-lock", GRAPH_ASSET_LOCK]
            if GRAPH_ASSETS_ENABLED
            else []
        ),
        exclude_sample_arguments=cli_repeated(
            "--exclude-truth-sample",
            config["pangenome"]["excluded_samples"],
        ),
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name build_pangenome_manifest \
          --job-key {PANGENOME_ID:q} \
          --run-id {RUN_ID:q} \
          --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} \
          --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} \
          --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} \
          --run-context {input.context:q} \
          --pangenome-manifest {output.pangenome_manifest:q} \
          --reference {input.reference:q} \
          --truth-profile {config[truth][primary]:q} \
          --snakemake-version {SNAKEMAKE_VERSION:q} \
          --execution-profile local \
          --random-seed {config[execution][random_seed]} \
          --threads 1 \
          --resource mem_mb=1024 \
          --param pangenome_id={PANGENOME_ID:q} \
          --param allele_namespace={config[pangenome][allele_namespace]:q} \
          --param excluded_samples={params.excluded_sample_csv:q} \
          --input {input.validated:q} \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.score_profile:q} \
          --input {input.reference:q} \
          --input {input.population:q} \
          --input {input.rule_source:q} \
          --input {input.rule_executor:q} \
          --input {input.script:q} \
          --input {input.environment:q} \
          --input {input.provenance_library:q} \
          {params.graph_input_arguments:q} \
          --output {output.panel:q} \
          --output {output.ledger:q} \
          --output {output.pangenome_manifest:q} \
          --upstream-manifest {input.validated_manifest:q} \
          --upstream-manifest {input.context_manifest:q} \
          {params.graph_upstream_arguments:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            --pangenome-id {PANGENOME_ID:q} \
            --backbone-id {config[reference][id]:q} \
            --reference {input.reference:q} \
            --population-source-id {config[pangenome][population_panel]:q} \
            --population-vcf {input.population:q} \
            --output-panel {output.panel:q} \
            --output-ledger {output.ledger:q} \
            --output-manifest {output.pangenome_manifest:q} \
            --allele-namespace {config[pangenome][allele_namespace]:q} \
            {params.exclude_sample_arguments:q} \
            {params.graph_script_arguments:q} \
          > {log:q} 2>&1
        """
