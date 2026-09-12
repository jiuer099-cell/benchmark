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


def optional_graph_asset(name):
    value = config["pangenome"]["graph_assets"].get(name)
    return [value] if value else []


if SHARED_ALIGNMENT_ENABLED:
    rule lock_shared_shortread_alignment:
        input:
            r1=config["sample"]["fastq_r1"],
            r2=config["sample"]["fastq_r2"],
            reference=config["reference"]["fasta"],
            bam=SHARED_ALIGNMENT_CONFIG["bam"],
            bai=SHARED_ALIGNMENT_CONFIG["bai"],
        output:
            lock=SHARED_ALIGNMENT_LOCK,
            rule_manifest=SHARED_ALIGNMENT_RULE_MANIFEST,
        log:
            f"{LOG_ROOT}/rules/lock_shared_shortread_alignment/{SAMPLE_ID}.log",
        benchmark:
            f"{BENCHMARK_ROOT}/rules/lock_shared_shortread_alignment/{SAMPLE_ID}.jsonl",
        conda:
            "../envs/core.yaml"
        shell:
            """
            {PYTHON_EXECUTABLE:q} {RULE_EXECUTOR:q} \\
              --rule-name lock_shared_shortread_alignment \\
              --job-key {SAMPLE_ID:q} --run-id {RUN_ID:q} \\
              --module-or-tool-id pgbench-core \\
              --snakefile-path workflow/rules/pangenome.smk \\
              --rule-source-path workflow/rules/pangenome.smk \\
              --script-or-wrapper-path workflow/scripts/validate_shared_alignment.py \\
              --config-snapshot {CONFIG_PATH:q} \\
              --score-profile {config[catalogs][score_weights]:q} \\
              --run-context {RESULTS_ROOT}/provenance/run-context.json \\
              --pangenome-manifest {RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml \\
              --reference {input.reference:q} --truth-profile {config[truth][primary]:q} \\
              --snakemake-version {SNAKEMAKE_VERSION:q} --execution-profile local \\
              --random-seed {config[execution][random_seed]} --threads 1 --resource mem_mb=1024 \\
              --input {input.r1:q} --input {input.r2:q} --input {input.reference:q} --input {input.bam:q} --input {input.bai:q} \\
              --output {output.lock:q} --manifest-output {output.rule_manifest:q} -- \\
              {PYTHON_EXECUTABLE:q} workflow/scripts/validate_shared_alignment.py \\
                --fastq-r1 {input.r1:q} --fastq-r2 {input.r2:q} --reference {input.reference:q} \\
                --bam {input.bam:q} --bai {input.bai:q} \\
                --source-fastq-sha256 {SHARED_ALIGNMENT_CONFIG[source_fastq_sha256]:q} \\
                --reference-sha256 {SHARED_ALIGNMENT_CONFIG[reference_sha256]:q} \\
                --aligner {SHARED_ALIGNMENT_CONFIG[aligner]:q} \\
                --aligner-version {SHARED_ALIGNMENT_CONFIG[aligner_version]:q} \\
                --command-sha256 {SHARED_ALIGNMENT_CONFIG[command_sha256]:q} \\
                --bam-sha256 {SHARED_ALIGNMENT_CONFIG[bam_sha256]:q} \\
                --bai-sha256 {SHARED_ALIGNMENT_CONFIG[bai_sha256]:q} \\
                --output {output.lock:q} > {log:q} 2>&1
            """


rule lock_graph_assets:
    input:
        validated=RESULTS_ROOT + "/provenance/config.validated.json",
        validated_manifest=VALIDATE_MANIFEST,
        context=RESULTS_ROOT + "/provenance/run-context.json",
        context_manifest=CONTEXT_MANIFEST,
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        reference=config["reference"]["fasta"],
        source_manifest=configured_graph_asset("manifest"),
        gbz=optional_graph_asset("gbz"),
        xg=optional_graph_asset("xg"),
        min_index=optional_graph_asset("min"),
        zipcodes=optional_graph_asset("zipcodes"),
        dist=optional_graph_asset("dist"),
        sample_list=optional_graph_asset("sample_list"),
        gfa=optional_graph_asset("gfa"),
        variation_calls=optional_graph_asset("variation_calls"),
        rule_source="workflow/rules/pangenome.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/validate_graph_assets.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    output:
        lock=GRAPH_ASSET_LOCK,
        rule_manifest=GRAPH_ASSET_RULE_MANIFEST,
    log:
        f"{LOG_ROOT}/rules/lock_graph_assets/{PANGENOME_ID}.log",
    benchmark:
        f"{BENCHMARK_ROOT}/rules/lock_graph_assets/{PANGENOME_ID}.jsonl",
    conda:
        "../envs/core.yaml"
    params:
        excluded_sample_csv=",".join(config["pangenome"]["excluded_samples"]),
        exclude_sample_arguments=cli_repeated(
            "--exclude-sample",
            config["pangenome"]["excluded_samples"],
        ),
        graph_input_arguments=cli_repeated(
            "--input",
            [
                config["pangenome"]["graph_assets"][name]
                for name in (
                    "manifest", "gbz", "xg", "min", "zipcodes", "dist", "sample_list", "gfa", "variation_calls"
                )
                if config["pangenome"]["graph_assets"].get(name)
            ],
        ),
        xg_script_arguments=(
            ["--xg", config["pangenome"]["graph_assets"]["xg"]]
            if config["pangenome"]["graph_assets"].get("xg")
            else []
        ),
        gbz_script_arguments=(
            ["--gbz", config["pangenome"]["graph_assets"]["gbz"]]
            if config["pangenome"]["graph_assets"].get("gbz")
            else []
        ),
        min_script_arguments=(
            ["--min", config["pangenome"]["graph_assets"]["min"]]
            if config["pangenome"]["graph_assets"].get("min")
            else []
        ),
        dist_script_arguments=(
            ["--dist", config["pangenome"]["graph_assets"]["dist"]]
            if config["pangenome"]["graph_assets"].get("dist")
            else []
        ),
        sample_list_script_arguments=(
            ["--sample-list", config["pangenome"]["graph_assets"]["sample_list"]]
            if config["pangenome"]["graph_assets"].get("sample_list")
            else []
        ),
        zipcodes_script_arguments=(
            ["--zipcodes", config["pangenome"]["graph_assets"]["zipcodes"]]
            if config["pangenome"]["graph_assets"].get("zipcodes")
            else []
        ),
        gfa_script_arguments=(
            ["--gfa", config["pangenome"]["graph_assets"]["gfa"]]
            if config["pangenome"]["graph_assets"].get("gfa")
            else []
        ),
        variation_calls_script_arguments=(
            [
                "--variation-calls",
                config["pangenome"]["graph_assets"]["variation_calls"],
            ]
            if config["pangenome"]["graph_assets"].get("variation_calls")
            else []
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
          --param graph_profile={GRAPH_PROFILE:q} \
          --param excluded_samples={params.excluded_sample_csv:q} \
          --input {input.validated:q} \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.score_profile:q} \
          --input {input.reference:q} \
          {params.graph_input_arguments:q} \
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
            --profile {GRAPH_PROFILE:q} \
            {params.gfa_script_arguments:q} \
            {params.variation_calls_script_arguments:q} \
            {params.gbz_script_arguments:q} \
            {params.xg_script_arguments:q} \
            {params.min_script_arguments:q} \
            {params.zipcodes_script_arguments:q} \
            {params.dist_script_arguments:q} \
            {params.sample_list_script_arguments:q} \
            --reference-path {config[pangenome][graph_assets][reference_path]:q} \
            {params.exclude_sample_arguments:q} \
            --output {output.lock:q} \
          > {log:q} 2>&1
        """


rule build_pangenome_manifest:
    input:
        validated=RESULTS_ROOT + "/provenance/config.validated.json",
        validated_manifest=VALIDATE_MANIFEST,
        context=RESULTS_ROOT + "/provenance/run-context.json",
        context_manifest=CONTEXT_MANIFEST,
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        reference=config["reference"]["fasta"],
        population=population_vcf_input,
        panel_provenance=config["catalogs"]["panel_provenance"],
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
        panel=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/panel.vcf",
        ledger=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/allele-ledger.tsv",
        pangenome_manifest=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml",
        rule_manifest=PANGENOME_RULE_MANIFEST,
    log:
        f"{LOG_ROOT}/rules/build_pangenome_manifest/{PANGENOME_ID}.log",
    benchmark:
        f"{BENCHMARK_ROOT}/rules/build_pangenome_manifest/{PANGENOME_ID}.jsonl",
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
          --input {input.panel_provenance:q} \
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
            --panel-provenance {input.panel_provenance:q} \
            --output-panel {output.panel:q} \
            --output-ledger {output.ledger:q} \
            --output-manifest {output.pangenome_manifest:q} \
            --allele-namespace {config[pangenome][allele_namespace]:q} \
            {params.exclude_sample_arguments:q} \
            {params.graph_script_arguments:q} \
          > {log:q} 2>&1
        """


if LEGACY_ADAPTER_PREPARATION_ENABLED:
    PANGENIE_CONTEXT = config["pangenome"].get("pangenie_private_context")
    if not isinstance(PANGENIE_CONTEXT, dict):
        raise WorkflowError(
            "PanGenie requires pangenome.pangenie_private_context in formal mode"
        )

    rule prepare_pangenie_private_panel:
        input:
            canonical_population=population_vcf_input,
            canonical_scoring_panel=(
                f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/"
                f"{SAMPLE_ID}.blinded.vcf"
            ),
            frozen_gbz=PANGENIE_CONTEXT["frozen_gbz"],
            frozen_sample_manifest=PANGENIE_CONTEXT["frozen_sample_manifest"],
            source_identity=PANGENIE_CONTEXT["source_identity"],
            source_graph=PANGENIE_CONTEXT["source_graph"],
            source_haplotype_manifest=PANGENIE_CONTEXT["source_haplotype_manifest"],
            source_phased_panel=PANGENIE_CONTEXT["source_phased_panel"],
            source_biallelic_panel=PANGENIE_CONTEXT["source_biallelic_panel"],
            source_biallelic_converter=PANGENIE_CONTEXT["source_biallelic_converter"],
            provenance=PANGENIE_CONTEXT["provenance"],
            pangenome_manifest=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml",
            pangenome_rule_manifest=PANGENOME_RULE_MANIFEST,
            challenge_rule_manifest=CHALLENGE_RULE_MANIFEST,
            context=RESULTS_ROOT + "/provenance/run-context.json",
            context_manifest=CONTEXT_MANIFEST,
            config=CONFIG_PATH,
            score_profile=config["catalogs"]["score_weights"],
            reference=config["reference"]["fasta"],
            rule_source="workflow/rules/pangenome.smk",
            rule_executor=RULE_EXECUTOR,
            script="workflow/scripts/prepare_pangenie_private_panel.py",
            environment=CORE_ENV_SPEC,
            provenance_library="workflow/scripts/pgbench_provenance.py",
        output:
            panel=PANGENIE_PRIVATE_PHASED_PANEL,
            projection=PANGENIE_PRIVATE_PROJECTION,
            gate=PANGENIE_PRIVATE_GATE,
            rule_manifest=PANGENIE_PRIVATE_RULE_MANIFEST,
        log:
            f"{LOG_ROOT}/rules/prepare_pangenie_private_panel/{SAMPLE_ID}.log",
        benchmark:
            f"{BENCHMARK_ROOT}/rules/prepare_pangenie_private_panel/{SAMPLE_ID}.jsonl",
        conda:
            "../envs/core.yaml"
        params:
            exclusions=cli_repeated(
                "--exclude-sample", config["pangenome"]["excluded_samples"]
            ),
        shell:
            """
            {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
              --rule-name prepare_pangenie_private_panel \
              --job-key {SAMPLE_ID:q} \
              --run-id {RUN_ID:q} \
              --module-or-tool-id pangenie \
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
              --resource mem_mb=4096 \
              --param source_cohort_id={PANGENIE_CONTEXT[source_cohort_id]:q} \
              --input {input.canonical_population:q} \
              --input {input.canonical_scoring_panel:q} \
              --input {input.frozen_gbz:q} \
              --input {input.frozen_sample_manifest:q} \
              --input {input.source_identity:q} \
              --input {input.source_graph:q} \
              --input {input.source_haplotype_manifest:q} \
              --input {input.source_phased_panel:q} \
              --input {input.source_biallelic_panel:q} \
              --input {input.source_biallelic_converter:q} \
              --input {input.provenance:q} \
              --input {input.pangenome_manifest:q} \
              --input {input.context:q} \
              --input {input.config:q} \
              --input {input.score_profile:q} \
              --input {input.reference:q} \
              --input {input.rule_source:q} \
              --input {input.rule_executor:q} \
              --input {input.script:q} \
              --input {input.environment:q} \
              --input {input.provenance_library:q} \
              --output {output.panel:q} \
              --output {output.projection:q} \
              --output {output.gate:q} \
              --upstream-manifest {input.pangenome_rule_manifest:q} \
              --upstream-manifest {input.challenge_rule_manifest:q} \
              --upstream-manifest {input.context_manifest:q} \
              --manifest-output {output.rule_manifest:q} \
              -- \
              {PYTHON_EXECUTABLE:q} {input.script:q} \
                --source-phased-panel {input.source_phased_panel:q} \
                --canonical-population-source {input.canonical_population:q} \
                --source-graph {input.source_graph:q} \
                --source-haplotype-manifest {input.source_haplotype_manifest:q} \
                --provenance {input.provenance:q} \
                --canonical-scoring-panel {input.canonical_scoring_panel:q} \
                --source-cohort-id {PANGENIE_CONTEXT[source_cohort_id]:q} \
                --frozen-gbz {input.frozen_gbz:q} \
                --frozen-sample-manifest {input.frozen_sample_manifest:q} \
                --source-identity {input.source_identity:q} \
                {params.exclusions:q} \
                --output-panel {output.panel:q} \
                --output-projection {output.projection:q} \
                --output-gate {output.gate:q} \
              > {log:q} 2>&1
            """
