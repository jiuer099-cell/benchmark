rule canonicalize_vcf:
    input:
        raw=external_raw_vcf,
        upstream=RESULTS_ROOT + "/provenance/rules/tool__{tool}__execute/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
        context=RESULTS_ROOT + "/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        pangenome_manifest=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_source="workflow/rules/normalization.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/normalize_sv_vcf.py",
        environment=CORE_ENV_SPEC,
        provenance_library=PROVENANCE_LIBRARY,
    output:
        vcf=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/calls.vcf",
        rule_manifest=RESULTS_ROOT + "/provenance/rules/canonicalize_vcf/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
    log: f"{LOG_ROOT}/rules/canonicalize_vcf/{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".log",
    benchmark: f"{BENCHMARK_ROOT}/rules/canonicalize_vcf/{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".jsonl",
    conda: "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name canonicalize_vcf --job-key {SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE} \
          --run-id {RUN_ID:q} --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} --run-context {input.context:q} \
          --pangenome-manifest {input.pangenome_manifest:q} --reference {input.reference:q} \
          --truth-profile {config[truth][primary]:q} --snakemake-version {SNAKEMAKE_VERSION:q} \
          --execution-profile local --random-seed {config[execution][random_seed]} \
          --threads 1 --resource mem_mb=1024 --wildcard tool={wildcards.tool:q} \
          --input {input.raw:q} --input {input.context:q} --input {input.config:q} \
          --input {input.score_profile:q} --input {input.pangenome_manifest:q} \
          --input {input.reference:q} --input {input.rule_source:q} \
          --input {input.rule_executor:q} --input {input.script:q} \
          --input {input.environment:q} --input {input.provenance_library:q} \
          --output {output.vcf:q} --upstream-manifest {input.upstream:q} \
          --manifest-output {output.rule_manifest:q} -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} --input-vcf {input.raw:q} \
            --reference {input.reference:q} --output-vcf {output.vcf:q} \
            --sample-name {SAMPLE_ID:q} > {log:q} 2>&1
        """


rule link_pangenome_alleles:
    input:
        canonical=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/calls.vcf",
        upstream=RESULTS_ROOT + "/provenance/rules/canonicalize_vcf/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
        ledger=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/allele-ledger.tsv",
        pangenome_manifest=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml",
        pangenome_upstream=PANGENOME_RULE_MANIFEST,
        context=RESULTS_ROOT + "/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        reference=config["reference"]["fasta"],
        rule_source="workflow/rules/normalization.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/link_pangenome_alleles.py",
        environment=CORE_ENV_SPEC,
        provenance_library=PROVENANCE_LIBRARY,
    output:
        vcf=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/linked.vcf",
        links=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/allele_links.tsv",
        rule_manifest=RESULTS_ROOT + "/provenance/rules/link_pangenome_alleles/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
    log: f"{LOG_ROOT}/rules/link_pangenome_alleles/{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".log",
    benchmark: f"{BENCHMARK_ROOT}/rules/link_pangenome_alleles/{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".jsonl",
    conda: "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name link_pangenome_alleles --job-key {SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE} \
          --run-id {RUN_ID:q} --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} --run-context {input.context:q} \
          --pangenome-manifest {input.pangenome_manifest:q} --reference {input.reference:q} \
          --truth-profile {config[truth][primary]:q} --snakemake-version {SNAKEMAKE_VERSION:q} \
          --execution-profile local --random-seed {config[execution][random_seed]} \
          --threads 1 --resource mem_mb=1024 --wildcard tool={wildcards.tool:q} \
          --input {input.canonical:q} --input {input.ledger:q} \
          --input {input.pangenome_manifest:q} --input {input.context:q} \
          --input {input.config:q} --input {input.score_profile:q} \
          --input {input.reference:q} --input {input.rule_source:q} \
          --input {input.rule_executor:q} --input {input.script:q} \
          --input {input.environment:q} --input {input.provenance_library:q} \
          --output {output.vcf:q} --output {output.links:q} \
          --upstream-manifest {input.upstream:q} --upstream-manifest {input.pangenome_upstream:q} \
          --manifest-output {output.rule_manifest:q} -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} --canonical-vcf {input.canonical:q} \
            --allele-ledger {input.ledger:q} --output-vcf {output.vcf:q} \
            --links-tsv {output.links:q} > {log:q} 2>&1
        """


rule materialize_all_sites:
    input:
        linked=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/linked.vcf",
        panel=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/{SAMPLE_ID}.blinded.vcf",
        upstream=RESULTS_ROOT + "/provenance/rules/link_pangenome_alleles/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
        context=RESULTS_ROOT + "/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        pangenome_manifest=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_source="workflow/rules/normalization.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/materialize_all_sites.py",
        environment=CORE_ENV_SPEC,
        provenance_library=PROVENANCE_LIBRARY,
    output:
        vcf=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/all-sites.vcf.gz",
        addressability=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/candidate-status.tsv",
        audit=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/addressability.audit.json",
        rule_manifest=RESULTS_ROOT + "/provenance/rules/materialize_all_sites/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
    log: f"{LOG_ROOT}/rules/materialize_all_sites/{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".log",
    benchmark: f"{BENCHMARK_ROOT}/rules/materialize_all_sites/{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".jsonl",
    conda: "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name materialize_all_sites --job-key {SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE} \
          --run-id {RUN_ID:q} --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} --run-context {input.context:q} \
          --pangenome-manifest {input.pangenome_manifest:q} --reference {input.reference:q} \
          --truth-profile {config[truth][primary]:q} --snakemake-version {SNAKEMAKE_VERSION:q} \
          --execution-profile local --random-seed {config[execution][random_seed]} \
          --threads 1 --resource mem_mb=1024 --wildcard tool={wildcards.tool:q} \
          --input {input.linked:q} --input {input.panel:q} --input {input.context:q} \
          --input {input.config:q} --input {input.score_profile:q} \
          --input {input.pangenome_manifest:q} --input {input.reference:q} \
          --input {input.rule_source:q} --input {input.rule_executor:q} \
          --input {input.script:q} --input {input.environment:q} \
          --input {input.provenance_library:q} --output {output.vcf:q} \
          --output {output.addressability:q} --output {output.audit:q} \
          --upstream-manifest {input.upstream:q} --manifest-output {output.rule_manifest:q} -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} --canonical-panel {input.panel:q} \
            --linked-query {input.linked:q} --output-vcf {output.vcf:q} \
            --addressability-tsv {output.addressability:q} --audit-json {output.audit:q} \
          > {log:q} 2>&1
        """


rule materialize_evaluation_query:
    input:
        panel=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/challenge/{SAMPLE_ID}.blinded.vcf",
        all_sites=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/all-sites.vcf.gz",
        candidate_status=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/candidate-status.tsv",
        upstream=RESULTS_ROOT + "/provenance/rules/materialize_all_sites/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
        context=RESULTS_ROOT + "/provenance/run-context.json",
        config=CONFIG_PATH,
        score_profile=config["catalogs"]["score_weights"],
        pangenome_manifest=f"{RESULTS_ROOT}/pangenome/{PANGENOME_ID}/manifest.yaml",
        reference=config["reference"]["fasta"],
        rule_source="workflow/rules/normalization.smk",
        rule_executor=RULE_EXECUTOR,
        script="workflow/scripts/materialize_evaluation_query.py",
        environment=CORE_ENV_SPEC,
        provenance_library=PROVENANCE_LIBRARY,
    output:
        vcf=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/evaluation-query.vcf.gz",
        audit=f"{RESULTS_ROOT}/{SAMPLE_ID}/{OFFICIAL_MODE}/" + "{tool}/canonical/evaluation-query.audit.json",
        rule_manifest=RESULTS_ROOT + "/provenance/rules/materialize_evaluation_query/" + f"{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".json",
    log: f"{LOG_ROOT}/rules/materialize_evaluation_query/{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".log",
    benchmark: f"{BENCHMARK_ROOT}/rules/materialize_evaluation_query/{SAMPLE_ID}." + "{tool}." + OFFICIAL_MODE + ".jsonl",
    conda: "../envs/core.yaml"
    shell:
        """
        {PYTHON_EXECUTABLE:q} {input.rule_executor:q} \
          --rule-name materialize_evaluation_query --job-key {SAMPLE_ID}.{wildcards.tool}.{OFFICIAL_MODE} \
          --run-id {RUN_ID:q} --module-or-tool-id pgbench-core \
          --snakefile-path {input.rule_source:q} --rule-source-path {input.rule_source:q} \
          --script-or-wrapper-path {input.script:q} --config-snapshot {input.config:q} \
          --score-profile {input.score_profile:q} --run-context {input.context:q} \
          --pangenome-manifest {input.pangenome_manifest:q} --reference {input.reference:q} \
          --truth-profile {config[truth][primary]:q} --snakemake-version {SNAKEMAKE_VERSION:q} \
          --execution-profile local --random-seed {config[execution][random_seed]} \
          --threads 1 --resource mem_mb=1024 --wildcard tool={wildcards.tool:q} \
          --input {input.panel:q} --input {input.all_sites:q} --input {input.candidate_status:q} \
          --input {input.context:q} --input {input.config:q} --input {input.score_profile:q} \
          --input {input.pangenome_manifest:q} --input {input.reference:q} \
          --input {input.rule_source:q} --input {input.rule_executor:q} \
          --input {input.script:q} --input {input.environment:q} \
          --input {input.provenance_library:q} --output {output.vcf:q} \
          --output {output.audit:q} --upstream-manifest {input.upstream:q} \
          --manifest-output {output.rule_manifest:q} -- \
          {PYTHON_EXECUTABLE:q} {input.script:q} --canonical-panel {input.panel:q} \
            --all-sites-vcf {input.all_sites:q} \
            --candidate-status-tsv {input.candidate_status:q} \
            --output-vcf {output.vcf:q} --audit-json {output.audit:q} \
          > {log:q} 2>&1
        """
