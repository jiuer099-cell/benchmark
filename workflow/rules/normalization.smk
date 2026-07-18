rule canonicalize_vcf:
    input:
        raw=external_raw_vcf,
        tool_rule_manifest=(
            "results/provenance/rules/tool__{tool}__execute/"
            f"{SAMPLE_ID}.{OFFICIAL_MODE}.json"
        ),
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_source="workflow/rules/normalization.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/normalize_sv_vcf.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    output:
        vcf=(
            f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            "{tool}/canonical/calls.vcf"
        ),
        rule_manifest=(
            "results/provenance/rules/canonicalize_vcf/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
    log:
        (
            f"logs/rules/canonicalize_vcf/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".log"
        ),
    benchmark:
        (
            f"benchmarks/rules/canonicalize_vcf/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".jsonl"
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        python {input.rule_executor:q} \
          --rule-name canonicalize_vcf \
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
          --resource mem_mb=1024 \
          --wildcard tool={wildcards.tool:q} \
          --param sample_name={SAMPLE_ID:q} \
          --input {input.raw:q} \
          --input {input.context:q} \
          --input {input.config:q} \
          --input {input.score_profile:q} \
          --input {input.pangenome_manifest:q} \
          --input {input.reference:q} \
          --input {input.rule_source:q} \
          --input {input.rule_executor:q} \
          --input {input.script:q} \
          --input {input.environment:q} \
          --input {input.provenance_library:q} \
          --output {output.vcf:q} \
          --upstream-manifest {input.tool_rule_manifest:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          python {input.script:q} \
            --input-vcf {input.raw:q} \
            --reference {input.reference:q} \
            --output-vcf {output.vcf:q} \
            --sample-name {SAMPLE_ID:q} \
          > {log:q} 2>&1
        """


rule link_pangenome_alleles:
    input:
        canonical=(
            f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            "{tool}/canonical/calls.vcf"
        ),
        canonical_rule_manifest=(
            "results/provenance/rules/canonicalize_vcf/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
        ledger=f"results/pangenome/{PANGENOME_ID}/allele-ledger.tsv",
        pangenome_manifest=f"results/pangenome/{PANGENOME_ID}/manifest.yaml",
        pangenome_rule_manifest=PANGENOME_RULE_MANIFEST,
        context="results/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        reference=config["reference"]["fasta"],
        rule_source="workflow/rules/normalization.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/link_pangenome_alleles.py",
        environment=CORE_ENV_SPEC,
        provenance_library="workflow/scripts/pgbench_provenance.py",
    output:
        vcf=(
            f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            "{tool}/canonical/linked.vcf"
        ),
        links=(
            f"results/{SAMPLE_ID}/{OFFICIAL_MODE}/"
            "{tool}/canonical/allele_links.tsv"
        ),
        rule_manifest=(
            "results/provenance/rules/link_pangenome_alleles/"
            f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json"
        ),
    log:
        (
            f"logs/rules/link_pangenome_alleles/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".log"
        ),
    benchmark:
        (
            f"benchmarks/rules/link_pangenome_alleles/{SAMPLE_ID}."
            "{tool}." + OFFICIAL_MODE + ".jsonl"
        ),
    conda:
        "../envs/core.yaml"
    shell:
        """
        python {input.rule_executor:q} \
          --rule-name link_pangenome_alleles \
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
          --resource mem_mb=1024 \
          --wildcard tool={wildcards.tool:q} \
          --input {input.canonical:q} \
          --input {input.ledger:q} \
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
          --output {output.vcf:q} \
          --output {output.links:q} \
          --upstream-manifest {input.canonical_rule_manifest:q} \
          --upstream-manifest {input.pangenome_rule_manifest:q} \
          --manifest-output {output.rule_manifest:q} \
          -- \
          python {input.script:q} \
            --canonical-vcf {input.canonical:q} \
            --allele-ledger {input.ledger:q} \
            --output-vcf {output.vcf:q} \
            --links-tsv {output.links:q} \
          > {log:q} 2>&1
        """
