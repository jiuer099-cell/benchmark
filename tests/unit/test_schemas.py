from __future__ import annotations

import csv
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")
jsonschema = pytest.importorskip("jsonschema")


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config"
SCHEMAS = ROOT / "workflow" / "schemas"


def _load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    assert isinstance(loaded, dict), f"{path} must contain a YAML mapping"
    return loaded


def _validator(schema_path: Path):
    schema = _load_yaml(schema_path)
    jsonschema.Draft202012Validator.check_schema(schema)
    return jsonschema.Draft202012Validator(
        schema,
        format_checker=jsonschema.FormatChecker(),
    )


def _assert_valid(schema_path: Path, instance: dict) -> None:
    errors = sorted(
        _validator(schema_path).iter_errors(instance),
        key=lambda error: list(error.absolute_path),
    )
    assert not errors, "\n".join(error.message for error in errors)


def _assert_invalid(schema_path: Path, instance: dict) -> None:
    assert list(_validator(schema_path).iter_errors(instance))


def test_all_registered_tool_manifests_match_schema() -> None:
    schema = SCHEMAS / "tool.schema.yaml"
    manifests = sorted((ROOT / "plugins").glob("*/tool.yaml"))
    assert manifests
    for path in manifests:
        _assert_valid(schema, _load_yaml(path))


def _valid_tool_manifest() -> dict:
    return {
        "interface_version": 1,
        "id": "my_genotyper",
        "name": "Synthetic external genotyper",
        "version": "0.1.0",
        "source": "external",
        "paradigm": "genotyping_only",
        "comparison_task": "panel_genotyping",
        "tasks": ["genotyping", "postprocess"],
        "inputs": {
            "reads": "bam",
            "reference": "required",
            "candidates": "required",
            "pangenome_manifest": "required",
            "pangenome_panel": "required",
            "graph_assets": {
                "required": False,
                "accepted_formats": ["vg", "xg_gbwt"],
            },
            "shared_alignment": {
                "required": True,
                "accepted_formats": ["bam"],
            },
            "read_sequences_auxiliary": False,
        },
        "capabilities": {
            "technology": ["pacbio_clr"],
            "svtypes": ["DEL", "INS", "DUP", "INV"],
            "phasing": False,
            "novel_discovery": False,
            "allele_linkage_output": "none",
            "allele_namespace": ["PGSV"],
        },
        "execution": {
            "runner": "plugins/my_genotyper/run.sh",
            "environment": "plugins/my_genotyper/envs/environment.yaml",
            "rule_registry": "plugins/my_genotyper/rule-registry.yaml",
            "trust_level": "untrusted",
            "sandbox_backend": "apptainer",
        },
        "supported_modes": {
            "caller_only_shared_alignment": {
                "required_inputs": [
                    "shared_alignment",
                    "reference",
                    "pangenome_manifest",
                    "pangenome_panel",
                    "candidate_panel",
                ],
                "optional_inputs": [],
                "forbidden_inputs": [
                    "canonical_fastq",
                    "original_input_bam",
                ],
                "billable_stages": ["genotype", "postprocess"],
            }
        },
        "outputs": {
            "vcf": "raw/calls.vcf.gz",
            "candidate_output_contract": "all_sites",
            "absence_semantics": "no_call",
        },
    }


def _valid_run_manifest() -> dict:
    digest = "a" * 64
    return {
        "manifest_schema_version": 1,
        "run_id": "synthetic_run",
        "sample_id": "HG002",
        "tool": {
            "id": "my_genotyper",
            "version": "0.1.0",
            "source": "external",
            "manifest_sha256": digest,
        },
        "command": ["plugins/my_genotyper/run.sh", "--threads", "4"],
        "parameters": {"threads": 4},
        "environment": {
            "conda_environment_export": "results/env/my_genotyper.yaml",
            "conda_lock_sha256": digest,
            "container_uri": None,
            "container_digest": None,
        },
        "inputs": [
            {
                "path": "prepared/HG002/all.fastq.gz",
                "size_bytes": 100,
                "sha256": digest,
                "read_only": True,
            }
        ],
        "reference": {"id": "grch38", "fasta_sha256": digest},
        "pangenome": {
            "manifest_id": "grch38_hprc_1kg_sv_v1",
            "manifest_sha256": digest,
            "panel_sha256": digest,
            "graph_assets": [],
            "allele_namespace": "PGSV",
        },
        "candidate_panel": {
            "id": "hg002_blinded_challenge_v1",
            "sha256": digest,
        },
        "official_score_mode": "end_to_end_from_reads",
        "alignment_kind": None,
        "alignment_group_id": None,
        "billable_stages": ["prepare_assets", "index", "map", "genotype"],
        "resources": {
            "threads": 4,
            "memory_mb": 8192,
            "disk_mb": 10240,
            "requested_runtime_minutes": 60,
        },
        "started_at": "2026-07-17T00:00:00Z",
        "finished_at": "2026-07-17T00:10:00Z",
        "exit_code": 0,
        "status": "completed",
        "outputs": [
            {
                "path": (
                    "results/HG002/end_to_end_from_reads/my_genotyper/raw/calls.vcf.gz"
                ),
                "size_bytes": 50,
                "sha256": digest,
                "read_only": False,
            }
        ],
    }


def _valid_pangenome_manifest() -> dict:
    digest = "b" * 64
    return {
        "schema_version": 1,
        "pangenome_id": "grch38_hprc_1kg_sv_v1",
        "backbone_reference": {
            "id": "grch38",
            "path": "resources/references/GRCh38_no_alt_plus_hs38d1_analysis_set.fasta",
            "sha256": digest,
        },
        "population_sources": [
            {
                "id": "hprc_1kg_grch38_sv",
                "path": "resources/pangenome/source.vcf.gz",
                "sha256": digest,
            }
        ],
        "allele_id_namespace": "PGSV",
        "truth_samples_excluded": ["HG002"],
        "panel_vcf": {
            "path": "resources/pangenome/panel.vcf.gz",
            "sha256": digest,
            "record_count": 12,
        },
        "allele_ledger": {
            "path": "resources/pangenome/allele-ledger.tsv.gz",
            "sha256": digest,
        },
        "graph_assets": None,
        "graph_build_recipe_sha256": None,
        "generated_at": "2026-07-17T00:00:00Z",
    }


def _valid_metric_document() -> dict:
    return {
        "schema_version": 1,
        "dictionary_id": "pgbench_metrics_v1",
        "tuple": {
            "run_id": "synthetic_run",
            "sample_id": "HG002",
            "tool_id": "my_genotyper",
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "giab_hg002_grch38_v5_0q",
            "score_profile": "pgbench_consensus_v3",
        },
        "records": [
            {
                "metric_id": "consensus.all_three_correct.count",
                "value_type": "count",
                "value": 75,
                "status": "defined",
                "evaluator": "fusion",
                "numerator": 75,
                "denominator": 100,
                "eligible_count": 100,
                "universe_id": "primary_truth_v5_0q",
                "undefined_reason": None,
                "parser_id": "parse_truvari_v1",
                "parser_source_field": "summary.f1",
                "provenance_manifest_id": "manifest-1",
                "strata": {},
                "confidence_interval": {
                    "lower": 0.70,
                    "upper": 0.80,
                    "level": 0.95,
                    "method": "stratified_variant_cluster_bootstrap",
                    "bootstrap_unit": "truth_variant_cluster",
                },
            }
        ],
    }


def test_all_json_schemas_are_draft_2020_12_and_well_formed() -> None:
    schema_paths = [
        CONFIG / "config.schema.yaml",
        *sorted(SCHEMAS.glob("*.schema.yaml")),
    ]
    assert {path.name for path in schema_paths} == {
        "config.schema.yaml",
        "metrics.schema.yaml",
        "pangenome-manifest.schema.yaml",
        "rule-registry.schema.yaml",
        "run-manifest.schema.yaml",
        "tool.schema.yaml",
    }
    for path in schema_paths:
        schema = _load_yaml(path)
        assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
        jsonschema.Draft202012Validator.check_schema(schema)


def test_example_config_validates_and_starts_without_unimplemented_tools() -> None:
    config = _load_yaml(CONFIG / "config.example.yaml")
    _assert_valid(CONFIG / "config.schema.yaml", config)
    assert config["tools"]["enabled"] == []
    assert config["tools"]["planned_tools"] == []
    assert config["score"]["produce_ranking"] is False
    assert config["score"]["final_score_name"] == "PGBenchConsensusScore"
    assert (
        config["execution"]["official_score_mode"]
        == "caller_only_shared_alignment"
    )
    assert config["development"]["synthetic_mode"] is False


def test_vg_end_to_end_example_config_matches_schema() -> None:
    config = _load_yaml(CONFIG / "config.vg.example.yaml")
    _assert_valid(CONFIG / "config.schema.yaml", config)
    assert config["execution"]["official_score_mode"] == "end_to_end_from_reads"
    assert config["sample"]["technology"] == "illumina_short_read"
    assert config["sample"]["fastq_r1"] is not None
    assert config["sample"]["fastq_r2"] is not None
    assert config["pangenome"]["graph_assets"]["profile"] == "vg_giraffe_shortread"
    assert config["external_plugins"] == [
        {"id": "vg", "manifest": "plugins/vg/tool.yaml"}
    ]


def test_pangenie_end_to_end_example_config_matches_schema() -> None:
    config = _load_yaml(CONFIG / "config.pangenie.example.yaml")
    _assert_valid(CONFIG / "config.schema.yaml", config)
    assert config["sample"]["technology"] == "illumina_short_read"
    assert config["sample"]["bam"] is None
    assert config["pangenome"]["build_graph_assets"] is False
    assert config["external_plugins"] == [
        {"id": "pangenie", "manifest": "plugins/pangenie/tool.yaml"}
    ]
    _assert_valid(
        SCHEMAS / "tool.schema.yaml",
        _load_yaml(ROOT / "plugins" / "pangenie" / "tool.yaml"),
    )


def test_caller_only_config_requires_non_null_bam() -> None:
    config = _load_yaml(CONFIG / "config.example.yaml")
    config["sample"]["bam"] = None
    _assert_invalid(CONFIG / "config.schema.yaml", config)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda config: config["tools"].update(enabled=["kanpig", "kanpig"]),
        lambda config: config["score"].update(produce_ranking=True),
        lambda config: config["reference"].pop("fasta"),
        lambda config: config["execution"].update(
            tracks=["caller_only_shared_alignment"],
            official_score_mode="end_to_end_from_reads",
        ),
        lambda config: config.pop("caller_only"),
    ],
)
def test_invalid_main_config_is_rejected(mutation) -> None:
    config = _load_yaml(CONFIG / "config.example.yaml")
    mutation(config)
    _assert_invalid(CONFIG / "config.schema.yaml", config)


def test_development_overrides_require_synthetic_mode() -> None:
    config = _load_yaml(CONFIG / "config.example.yaml")
    config["development"]["canonical_fastq"] = "tests/fixtures/all.fastq.gz"
    config["development"]["population_vcf"] = "tests/fixtures/population.vcf"
    config["development"]["truth_vcf"] = "tests/fixtures/truth.vcf"
    config["development"]["benchmark_bed"] = "tests/fixtures/benchmark.bed"
    config["development"]["allow_missing_bam"] = True
    _assert_invalid(CONFIG / "config.schema.yaml", config)

    config["development"]["synthetic_mode"] = True
    _assert_valid(CONFIG / "config.schema.yaml", config)


def test_sample_table_has_the_canonical_hg002_bam() -> None:
    with (CONFIG / "samples.example.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["sample_id"] == "HG002"
    assert rows[0]["technology"] == "pacbio_clr"
    assert rows[0]["bam"].endswith(
        "HG002.GRCh38.bam"
    )


def test_truth_and_stratification_profiles_are_frozen() -> None:
    truthsets = _load_yaml(CONFIG / "truthsets.yaml")
    primary = truthsets["truthsets"]["giab_hg002_grch38_v5_0q"]
    assert primary["role"] == "primary"
    assert primary["truth_status"] == "draft"
    assert primary["legacy"] is False
    assert primary["source_page"].endswith(
        "/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/"
    )
    assert truthsets["reference_id"] == "grch38"
    assert {
        name: asset["sha256"]
        for name, asset in primary["files"].items()
    } == {
        "vcf": "d66b2d2496ff5418763813d3195599007dbff950514d10e5bc27b6c8b76e34b8",
        "vcf_index": (
            "57898d17e44ecd5c286aeb863f27e613244396a3ad4b04ca79c8add647b8a81d"
        ),
        "benchmark_bed": (
            "2f75ce942e1dd9e1a443e4e04fac640104aae5e23fba24d7d3f2c1f1f9985b00"
        ),
    }
    legacy = truthsets["truthsets"][
        "giab_hg002_grch38_t2tq100_v0_9_sv_legacy"
    ]
    assert legacy["role"] == "legacy"
    assert legacy["legacy"] is True
    assert legacy["superseded_by"] == "giab_hg002_grch38_v5_0q"

    stratifications = _load_yaml(CONFIG / "stratifications.yaml")
    required = {
        name
        for name, item in stratifications["strata"].items()
        if item["required_for_score"]
    }
    assert {
        "difficult",
        "segmental_duplication",
        "tandem_repeat",
        "homopolymer",
        "low_complexity",
        "low_mappability",
        "gc_extreme",
        "coding_exonic",
        "medically_relevant",
    } <= required


def test_external_genotyper_manifest_validates() -> None:
    _assert_valid(SCHEMAS / "tool.schema.yaml", _valid_tool_manifest())


def test_genotyping_only_tool_cannot_submit_variant_sites_contract() -> None:
    manifest = _valid_tool_manifest()
    manifest["outputs"]["candidate_output_contract"] = "variant_sites"
    manifest["outputs"]["absence_semantics"] = "hom_ref"
    _assert_invalid(SCHEMAS / "tool.schema.yaml", manifest)


def test_external_tool_accepts_runner_or_workflow_module() -> None:
    manifest = _valid_tool_manifest()
    manifest["execution"].pop("runner")
    manifest["execution"]["module"] = "plugins/my_genotyper/Snakefile"
    manifest["execution"]["trust_level"] = "trusted_reviewed"
    _assert_valid(SCHEMAS / "tool.schema.yaml", manifest)

    manifest["execution"].pop("module")
    _assert_invalid(SCHEMAS / "tool.schema.yaml", manifest)


def test_external_tool_rejects_runner_and_module_together() -> None:
    manifest = _valid_tool_manifest()
    manifest["execution"]["module"] = "plugins/my_genotyper/Snakefile"
    _assert_invalid(SCHEMAS / "tool.schema.yaml", manifest)


def test_external_tool_rejects_output_path_traversal() -> None:
    manifest = _valid_tool_manifest()
    manifest["outputs"]["vcf"] = "../precomputed.vcf"
    _assert_invalid(SCHEMAS / "tool.schema.yaml", manifest)


def test_rule_registry_accepts_normative_template_ids() -> None:
    registry = {
        "schema_version": 1,
        "rules": [
            {
                "id": "tool__{tool}__genotype",
                "description": "Run a registered genotyper.",
                "phase": "tool_execution",
                "kind": "tool_template",
                "inputs": ["shared_alignment", "candidate_panel"],
                "outputs": ["raw_vcf"],
                "upstream": ["tool__{tool}__preflight"],
                "audit_scope": "pre_score",
                "manifest_required": True,
                "log_required": True,
                "benchmark_required": True,
                "implementation_status": "partial",
                "source_path": "workflow/rules/tools.smk",
                "environment": None,
                "resource_class": "tool_default",
                "active_when": {"task": "genotyping"},
            }
        ],
    }
    _assert_valid(SCHEMAS / "rule-registry.schema.yaml", registry)

    registry["rules"][0]["id"] = "tool/genotype"
    _assert_invalid(SCHEMAS / "rule-registry.schema.yaml", registry)


def test_run_manifest_validates_and_rejects_bad_hashes() -> None:
    manifest = _valid_run_manifest()
    _assert_valid(SCHEMAS / "run-manifest.schema.yaml", manifest)
    manifest["outputs"][0]["sha256"] = "not-a-sha256"
    _assert_invalid(SCHEMAS / "run-manifest.schema.yaml", manifest)


def test_completed_run_requires_zero_exit_and_output() -> None:
    manifest = _valid_run_manifest()
    manifest["exit_code"] = 1
    manifest["outputs"] = []
    _assert_invalid(SCHEMAS / "run-manifest.schema.yaml", manifest)


def test_pangenome_manifest_requires_hg002_truth_exclusion() -> None:
    manifest = _valid_pangenome_manifest()
    _assert_valid(SCHEMAS / "pangenome-manifest.schema.yaml", manifest)
    manifest["truth_samples_excluded"] = ["HG001"]
    _assert_invalid(SCHEMAS / "pangenome-manifest.schema.yaml", manifest)


def test_pangenome_manifest_accepts_content_locked_graph_assets() -> None:
    manifest = _valid_pangenome_manifest()
    digest = "c" * 64
    asset = {
        "path": "resources/pangenome/graph/graph.gbz",
        "sha256": digest,
        "size_bytes": 100,
    }
    manifest["graph_assets"] = {
        "profile": "vg_legacy_xg",
        "manifest": {
            **asset,
            "path": "resources/pangenome/graph/graph-assets.lock.yaml",
        },
        "lock_manifest": {
            **asset,
            "path": "results/pangenome/graph-assets.lock.yaml",
        },
        "asset_root": "resources/pangenome/graph",
        "reference_path": "GRCh38",
        "excluded_samples": ["HG002", "NA24385"],
        "gbz": asset,
        "xg": {**asset, "path": "resources/pangenome/graph/graph.xg"},
        "min": {**asset, "path": "resources/pangenome/graph/graph.min"},
        "dist": {**asset, "path": "resources/pangenome/graph/graph.dist"},
        "sample_list": {
            **asset,
            "path": "resources/pangenome/graph/samples.txt",
        },
        "sample_count": 2,
    }
    _assert_valid(SCHEMAS / "pangenome-manifest.schema.yaml", manifest)

    manifest["graph_assets"]["gbz"]["sha256"] = "unlocked"
    _assert_invalid(SCHEMAS / "pangenome-manifest.schema.yaml", manifest)


def test_metrics_document_validates() -> None:
    _assert_valid(SCHEMAS / "metrics.schema.yaml", _valid_metric_document())


def test_metrics_reject_unknown_ids_and_negative_counts() -> None:
    document = _valid_metric_document()
    document["records"][0]["metric_id"] = "temporary.parser.metric"
    _assert_invalid(SCHEMAS / "metrics.schema.yaml", document)

    document = _valid_metric_document()
    document["records"][0]["value"] = -1
    _assert_invalid(SCHEMAS / "metrics.schema.yaml", document)


def test_undefined_metric_requires_null_value_and_reason() -> None:
    document = _valid_metric_document()
    record = document["records"][0]
    record["status"] = "undefined"
    record["value"] = None
    record["undefined_reason"] = "No truth-positive records in frozen stratum."
    record["confidence_interval"] = None
    _assert_valid(SCHEMAS / "metrics.schema.yaml", document)

    record["undefined_reason"] = None
    _assert_invalid(SCHEMAS / "metrics.schema.yaml", document)


def test_consensus_profile_has_three_equal_votes_and_no_weights() -> None:
    profile = _load_yaml(CONFIG / "consensus_scoring.yaml")
    assert profile["profile"]["id"] == "pgbench_consensus_v3"
    assert profile["consensus"]["evaluators"] == ["truvari", "aardvark", "vcfdist"]
    assert "layers" not in profile
    assert "evaluator_weights" not in profile["consensus"]


def test_score_metrics_match_dictionary_and_schema_enum() -> None:
    dictionary = _load_yaml(CONFIG / "metric_dictionary.yaml")
    dictionary_ids = set(dictionary["metrics"])
    metric_schema = _load_yaml(SCHEMAS / "metrics.schema.yaml")
    schema_ids = set(metric_schema["$defs"]["metricId"]["enum"])
    assert dictionary_ids == schema_ids

    assert {
        "consensus.all_three_correct.count",
        "consensus.exactly_two_correct.count",
        "consensus.exactly_one_correct.count",
        "consensus.none_correct.count",
    } <= dictionary_ids


def test_no_weight_profile_remains() -> None:
    assert not (CONFIG / "score_weights.yaml").exists()
