def population_vcf_input(wildcards):
    path = config.get("development", {}).get("population_vcf")
    if path:
        return path
    raise WorkflowError(
        "production population panel resolution is not implemented yet; "
        "configure a frozen catalog asset in Phase 2"
    )


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
          --output {output.panel:q} \
          --output {output.ledger:q} \
          --output {output.pangenome_manifest:q} \
          --upstream-manifest {input.validated_manifest:q} \
          --upstream-manifest {input.context_manifest:q} \
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
          > {log:q} 2>&1
        """
