def truth_vcf_input(wildcards):
    path = config.get("development", {}).get("truth_vcf")
    if path:
        return path
    path = config.get("evaluation", {}).get("truth_vcf")
    if path:
        return path
    raise WorkflowError("evaluation.truth_vcf is required in formal mode")


rule build_blinded_challenge_panel:
    input:
        panel=f"results/pangenome/{PANGENOME_ID}/panel.vcf",
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        pangenome_rule_manifest=PANGENOME_RULE_MANIFEST,
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        reference=config["reference"]["fasta"],
        truth=truth_vcf_input,
        rule_source="workflow/rules/panel.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/build_challenge_panel.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    output:
        vcf=(
            f"results/pangenome/{PANGENOME_ID}/challenge/"
            f"{SAMPLE_ID}.blinded.vcf"
        ),
        ledger=(
            f"results/pangenome/{PANGENOME_ID}/challenge/"
            f"{SAMPLE_ID}.hidden.tsv"
        ),
        audit=(
            f"results/pangenome/{PANGENOME_ID}/challenge/"
            f"{SAMPLE_ID}.audit.json"
        ),
        rule_manifest=CHALLENGE_RULE_MANIFEST,
    log:
        f"logs/rules/build_blinded_challenge_panel/{SAMPLE_ID}.log",
    benchmark:
        f"benchmarks/rules/build_blinded_challenge_panel/{SAMPLE_ID}.jsonl",
    conda:
        "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name build_blinded_challenge_panel \
          --job-key {SAMPLE_ID:q} \
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
          --resource mem_mb=1024 \
          --param seed={config[execution][random_seed]} \
          --input {input.panel:q} \
          --input {input.pangenome_manifest:q} \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.score_profile:q} \
          --input {input.reference:q} \
          --input {input.truth:q} \
          --input {input.rule_source:q} \
          --input {input.rule_executor:q} \
          --input {input.script:q} \
          --input {input.environment:q} \
          --input {input.provenance_library:q} \
          --output {output.vcf:q} \
          --output {output.ledger:q} \
          --output {output.audit:q} \
          --upstream-manifest {input.pangenome_rule_manifest:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} \
            --panel-vcf {input.panel:q} \
            --truth-vcf {input.truth:q} \
            --output-vcf {output.vcf:q} \
            --hidden-ledger {output.ledger:q} \
            --audit-json {output.audit:q} \
            --seed {config[execution][random_seed]} \
          > {log:q} 2>&1
        """
