REPORT_RULE_MANIFESTS = [
    semantic_rule_manifest(
        "render_report",
        f"{SAMPLE_ID}.{settings['tool_id']}.{OFFICIAL_MODE}",
    )
    for settings in EXTERNAL_SETTINGS
]
FINALIZER_RULE_MANIFESTS = [
    semantic_rule_manifest(
        "finalize_score_provenance",
        f"{SAMPLE_ID}.{settings['tool_id']}.{OFFICIAL_MODE}",
    )
    for settings in EXTERNAL_SETTINGS
]


rule render_report:
    input:
        score="results/summary/{tool}/score.json",
        package="results/summary/{tool}/score-package.json",
        finalizer_rule_manifest=(
            "results/provenance/rules/finalize_score_provenance/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_executor=RULE_EXECUTOR,
        rule_source="workflow/rules/report.smk",
        core_env="workflow/envs/core.yaml",
        script="workflow/scripts/render_report.py",
        provenance_library=PROVENANCE_LIBRARY,
    output:
        html="results/report/tool_cards/{tool}.html",
        score_tsv="results/summary/{tool}/score.tsv",
        points="results/summary/{tool}/point_breakdown.tsv",
        rule_manifest=(
            "results/provenance/rules/render_report/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
    log:
        (
            f"logs/rules/render_report/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".log"
        ),
    benchmark:
        (
            f"benchmarks/rules/render_report/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".jsonl"
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name render_report \
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
          --input {input.score:q} \
          --input {input.package:q} \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.score_profile:q} \
          --input {input.pangenome_manifest:q} \
          --input {input.reference:q} \
          --input {input.rule_executor:q} \
          --input {input.rule_source:q} \
          --input {input.core_env:q} \
          --input {input.script:q} \
          --input {input.provenance_library:q} \
          --output {output.html:q} \
          --output {output.score_tsv:q} \
          --output {output.points:q} \
          --upstream-manifest {input.finalizer_rule_manifest:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            --score-json {input.score:q} \
            --score-package {input.package:q} \
            --finalizer-manifest {input.finalizer_rule_manifest:q} \
            --html {output.html:q} \
            --score-tsv {output.score_tsv:q} \
            --point-breakdown-tsv {output.points:q} \
          > {log:q} 2>&1
        """


rule render_report_index:
    input:
        cards=REPORT_CARDS,
        report_rule_manifests=REPORT_RULE_MANIFESTS,
        scores=SCORE_JSONS,
        metrics=METRICS_JSONS,
        packages=FINAL_SCORE_PACKAGES,
        tool_manifests=TOOL_MANIFESTS,
        finalizer_rule_manifests=FINALIZER_RULE_MANIFESTS,
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_executor=RULE_EXECUTOR,
        rule_source="workflow/rules/report.smk",
        core_env="workflow/envs/core.yaml",
        script="workflow/scripts/render_report_index.py",
        aggregate_script="workflow/scripts/aggregate_run_summary.py",
        provenance_library=PROVENANCE_LIBRARY,
    output:
        html="results/report/index.html",
        score_tsv="results/summary/score.tsv",
        points_tsv="results/summary/point_breakdown.tsv",
        metrics_tsv="results/summary/metrics.long.tsv",
        metrics_json="results/summary/metrics.json",
        rule_manifest=(
            "results/provenance/rules/render_report_index/index.json"
        ),
    log:
        "logs/rules/render_report_index/index.log",
    benchmark:
        "benchmarks/rules/render_report_index/index.jsonl",
    params:
        card_args=cli_repeated("--card", REPORT_CARDS),
        score_args=cli_repeated("--score-json", SCORE_JSONS),
        metrics_args=cli_repeated("--metrics-json", METRICS_JSONS),
        package_args=cli_repeated("--score-package", FINAL_SCORE_PACKAGES),
        finalizer_args=cli_repeated(
            "--finalizer-manifest", FINALIZER_RULE_MANIFESTS
        ),
        tool_manifest_args=cli_repeated(
            "--tool-manifest", TOOL_MANIFESTS
        ),
        provenance_input_args=cli_repeated(
            "--input",
            [
                *REPORT_CARDS,
                *REPORT_RULE_MANIFESTS,
                *SCORE_JSONS,
                *METRICS_JSONS,
                *FINAL_SCORE_PACKAGES,
                *TOOL_MANIFESTS,
                *FINALIZER_RULE_MANIFESTS,
                "results/provenance/run-context.json",
                CONFIG_PATH,
                config["catalogs"]["score_weights"],
                f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
                config["reference"]["fasta"],
                RULE_EXECUTOR,
                "workflow/rules/report.smk",
                "workflow/envs/core.yaml",
                "workflow/scripts/render_report_index.py",
                "workflow/scripts/aggregate_run_summary.py",
                PROVENANCE_LIBRARY,
            ],
        ),
        upstream_args=cli_repeated(
            "--upstream-manifest",
            [*REPORT_RULE_MANIFESTS, *FINALIZER_RULE_MANIFESTS],
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name render_report_index \
          --job-key index \
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
          {params.provenance_input_args:q} \
          --output {output.html:q} \
          --output {output.score_tsv:q} \
          --output {output.points_tsv:q} \
          --output {output.metrics_tsv:q} \
          --output {output.metrics_json:q} \
          {params.upstream_args:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            {params.card_args:q} \
            {params.score_args:q} \
            {params.metrics_args:q} \
            {params.package_args:q} \
            {params.finalizer_args:q} \
            {params.tool_manifest_args:q} \
            --output {output.html:q} \
            --output-score-tsv {output.score_tsv:q} \
            --output-point-tsv {output.points_tsv:q} \
            --output-metrics-tsv {output.metrics_tsv:q} \
            --output-metrics-json {output.metrics_json:q} \
          > {log:q} 2>&1
        """
