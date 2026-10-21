"""Performance-only regression checks against the frozen student build."""

import ast
import copy
import csv
import gzip
import hashlib
import itertools
import json
import os
import random
import subprocess
import sys
import types
from argparse import Namespace
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import run_viharmony as harmony
import test_viharmony as existing_tests
from reporting_test_paths import artifact, artifacts

write_tsv = existing_tests.write_tsv

REFERENCE_COMMIT = "d4b90f3c235eae41b03a24ea4fbfd62b5febdb1c"


def exhaustive_join(evidence_rows, final_to_parent, regions, refined):
    """Original all-pairs join, deliberately kept unoptimized as a test oracle."""
    known_final_ids = set(final_to_parent)
    result = defaultdict(list)
    for row in evidence_rows:
        exact = row.get("sequence_id", "").strip()
        parent = row.get("parent_sequence_id", "").strip() or exact
        record_type = row.get("record_type", "").strip().lower()
        interval = harmony.parse_optional_interval(row.get("coordinates", ""))
        for final_id, final_parent in final_to_parent.items():
            if harmony.evidence_applies_to_final(
                exact, parent, record_type, interval, final_id, final_parent,
                regions[final_id].get("record_type", "").strip(),
                harmony.parse_optional_interval(regions[final_id].get("coordinates", "")),
                len(refined[final_id]), known_final_ids,
            ):
                result[final_id].append(row)
    return result


def assert_same_join(rows, parents, regions, refined):
    expected = exhaustive_join(rows, parents, regions, refined)
    actual = harmony.index_evidence_by_final(rows, parents, regions, refined)
    assert list(actual) == list(expected)
    # Include row identity and multiplicity, not just value equality.
    for final_id in expected:
        assert [id(row) for row in actual[final_id]] == [id(row) for row in expected[final_id]]


def test_all_applicability_branches_preserve_row_identity_and_order():
    parents = {"a_child": "a", "b_child": "b", "a_second": "a", "intact": "intact",
               "a_full": "a", "intact_child": "intact", "b_full": "b", "blank": ""}
    regions = {
        final_id: {"record_type": "provirus", "coordinates": "11-50"}
        for final_id in parents
    }
    regions["a_second"] = {"record_type": "viral_region", "coordinates": "51-90"}
    for final_id in ("intact", "a_full", "b_full", "blank"):
        regions[final_id] = {"record_type": "input_contig", "coordinates": ""}
    refined = {final_id: "A" * 100 for final_id in parents}
    rows = [
        dict(sequence_id=exact, parent_sequence_id=parent, record_type=record_type, coordinates=coords)
        for exact, parent, record_type, coords in itertools.product(
            ("a", "b", "intact", "a_child", "missing", " a ", ""),
            ("a", "b", "intact", "missing", " a ", "", " "),
            ("input_contig", "provirus", "viral_region", "short_hallmark_region", " PROVIRUS ", ""),
            ("", "NA", "bad", "1-100", "11-50", "50-51", "91-100", "0-10"),
        )
    ]
    rows.extend([rows[0], dict(rows[0]), {}])
    assert_same_join(rows, parents, regions, refined)


@pytest.mark.parametrize("seed", range(10))
def test_randomized_join_matches_all_pairs(seed):
    rng = random.Random(seed)
    parents = {f"f{i}": rng.choice(["p0", "p1", "p2", "f0", ""]) for i in range(25)}
    regions = {key: {"record_type": rng.choice(["input_contig", "provirus", "viral_region"]),
                     "coordinates": rng.choice(["", "1-50", "30-100"])} for key in parents}
    refined = {key: "A" * rng.choice([50, 100]) for key in parents}
    ids = [*parents, "p0", "p1", "p2", "absent", ""]
    rows = [dict(sequence_id=rng.choice(ids), parent_sequence_id=rng.choice(ids),
                 record_type=rng.choice(["", "provirus", "viral_region", "short_hallmark_region"]),
                 coordinates=rng.choice(["", "1-50", "1-100", "30-100"])) for _ in range(200)]
    assert_same_join(rows, parents, regions, refined)


def test_empty_join():
    assert_same_join([], {}, {}, {})
    assert_same_join([dict(sequence_id="missing")], {}, {}, {})


def test_index_avoids_unrelated_final_sequences():
    count = 2000
    parents = {f"f{i}": f"p{i}" for i in range(count)}
    regions = {key: {"record_type": "provirus", "coordinates": "1-50"} for key in parents}
    refined = {key: "A" * 50 for key in parents}
    rows = [dict(sequence_id=parent) for parent in parents.values()]
    with patch.object(harmony, "evidence_applies_to_final", wraps=harmony.evidence_applies_to_final) as predicate:
        result = harmony.index_evidence_by_final(rows, parents, regions, refined)
    assert len(result) == count
    assert predicate.call_count == count  # Original algorithm: count * count.


@pytest.fixture(scope="module")
def frozen_reference():
    try:
        source = subprocess.check_output(
            ["git", "show", f"{REFERENCE_COMMIT}:bin/run_viharmony.py"], cwd=ROOT,
            stderr=subprocess.PIPE, text=True, encoding="utf-8",
        )
    except (OSError, subprocess.CalledProcessError):
        pytest.skip(f"Frozen-build replay requires Git object {REFERENCE_COMMIT}")
    module = types.ModuleType("_viharmony_frozen_student_reference")
    sys.modules[module.__name__] = module
    exec(compile(source, "<frozen student viHARMONY>", "exec"), module.__dict__)
    yield module, source
    del sys.modules[module.__name__]


def compare_outputs(reference_prefix, actual_prefix):
    """Keep the frozen oracle, with explicit expectations for the one TE fixture.

    New report files/columns are checked separately. All original metadata fields
    remain exact except the explicitly named retro fixture's approved outcomes.
    No arbitrary TE-row or origin-field exclusion is allowed.
    """
    def table(path):
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return reader.fieldnames, list(reader)

    def sequences(path):
        result = {}; sid = None
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                sid = line[1:]
                assert sid not in result
                result[sid] = ""
            elif line:
                result[sid] += line
        return result

    expected = {p.name: p for p in reference_prefix.parent.glob(reference_prefix.name + ".*")}
    actual = artifacts(actual_prefix)
    new_suffixes = ("taxonomy_reporting_audit.tsv.gz", "candidate_routing_audit.tsv.gz",
        "TE_candidates.fasta", "TE_candidates.summary.tsv", "viral_candidates.fasta", "viral_candidates.summary.tsv")
    actual_manifest = json.loads(actual[actual_prefix.name + ".harmonizer_manifest.json"].read_text())
    report_labels = {"taxonomy", "checkv_report", "vcontact3_report", "review_report"}
    report_labels.update(k for k in actual_manifest["outputs"] if k.startswith("review_reason_"))
    report_names = {Path(actual_manifest["outputs"][k]["file"]).name for k in report_labels}
    removed_name = reference_prefix.name + ".review_candidates.fasta"
    assert set(actual) == (set(expected) - {removed_name}) | {actual_prefix.name + "." + suffix for suffix in new_suffixes} | report_names
    reference_manifest = json.loads(expected[reference_prefix.name + ".harmonizer_manifest.json"].read_text())
    rna = reference_manifest["input_type"] == "rna"
    columns, baseline = table(expected[reference_prefix.name + ".final_metadata.tsv"])
    actual_columns, native = table(actual[actual_prefix.name + ".final_metadata.tsv"])
    assert actual_columns == columns + harmony.te_routing.METADATA_COLUMNS
    adjusted = []
    for original in baseline:
        row = dict(original)
        if row["final_sequence_id"] == "retro":
            # This fixture has only generic TE annotation (no complete direct
            # architecture), plus qualified Deep6/viCAT. It gains eligibility.
            assert row["viral_decision"] == "ambiguous_review"
            row["viral_decision"] = "provisional_viral" if rna else "retained_viral"
            row["viral_entity_interpretation"] = "provisional_viral_contig" if rna else "viral_contig"
            if rna:
                row["viral_confidence"] = "provisional"  # Original RNA floor still applies.
        adjusted.append(row)
    assert [{k:r[k] for k in columns} for r in native] == adjusted
    for row,original in zip(native,baseline,strict=True):
        assert row["pre_TE_routing_decision"] == original["viral_decision"]
        if row["final_sequence_id"] != "retro":
            assert row["viral_decision"] == original["viral_decision"]

    primary = {r["final_sequence_id"] for r in adjusted if r["viral_decision"] in harmony.PRIMARY_VIRAL_DECISIONS}
    provisional = {r["final_sequence_id"] for r in adjusted if r["viral_decision"] in harmony.PROVISIONAL_VIRAL_DECISIONS}
    all_sequences = sequences(expected[reference_prefix.name + ".all_candidates.fasta"])
    _, mapping = table(expected[reference_prefix.name + ".sequence_map.tsv"])
    original_ids = {r["final_sequence_id"]:r["final_original_id"] for r in mapping}
    for suffix,ids,original in (
        ("final.normalized.fasta",primary,False),("final.original_ids.fasta",primary,True),
        ("database_candidates.fasta",primary,False),
        ("provisional.normalized.fasta",provisional,False),("provisional.original_ids.fasta",provisional,True),
        ("viral_candidates.fasta",primary|provisional,False),("TE_candidates.fasta",set(),False)):
        want={(original_ids[sid] if original else sid):seq for sid,seq in all_sequences.items() if sid in ids}
        observed=sequences(actual[actual_prefix.name + "." + suffix])
        assert list(observed.items()) == list(want.items()),suffix

    for suffix,ids in (("database_candidates.tsv",primary),("provisional_metadata.tsv",provisional)):
        header,observed=table(actual[actual_prefix.name + "." + suffix])
        assert header == actual_columns
        assert [{k:r[k] for k in columns} for r in observed] == [r for r in adjusted if r["final_sequence_id"] in ids]
    review_columns,old_review=table(expected[reference_prefix.name + ".review_queue.tsv"])
    _,new_review=table(actual[actual_prefix.name + ".review_queue.tsv"])
    adjusted_index={r["final_sequence_id"]:r for r in adjusted}
    wanted_review=[]
    for row in old_review:
        result=dict(row,**adjusted_index[row["final_sequence_id"]])
        if row["final_sequence_id"] == "retro" and rna:
            result["review_reasons"]="single_qualified_viral_tool,"+row["review_reasons"]
        wanted_review.append(result)
    assert [{k:r[k] for k in review_columns} for r in new_review] == wanted_review
    header,old_disposition=table(expected[reference_prefix.name + ".sequence_disposition.tsv"])
    _,new_disposition=table(actual[actual_prefix.name + ".sequence_disposition.tsv"])
    for row in old_disposition:
        if row["normalized_name"] == "retro":
            row["disposition"]="provisional_viral_after_harmony" if rna else "retained_contig"
            row["final_sequence_ids"]="NA" if rna else "retro"
            row["provisional_sequence_ids"]="retro" if rna else "NA"
    assert new_disposition == old_disposition

    handled=("final_metadata.tsv","database_candidates.tsv","provisional_metadata.tsv", "review_queue.tsv",
        "sequence_disposition.tsv","final.normalized.fasta","final.original_ids.fasta","database_candidates.fasta",
        "provisional.normalized.fasta","provisional.original_ids.fasta","harmonizer_manifest.json", "review_candidates.fasta")
    for name,left in expected.items():
        if name in {reference_prefix.name + "." + suffix for suffix in handled}:
            continue
        right=actual[name]
        a,b=left.read_bytes(),right.read_bytes()
        if name.endswith(".gz"):
            a,b=gzip.decompress(a),gzip.decompress(b)
        assert a==b,name
    native_manifest=json.loads(actual[actual_prefix.name + ".harmonizer_manifest.json"].read_text())
    # Verify every output checksum, including the six additions, before removing
    # only approved manifest fields from this historical-build comparison.
    for item in native_manifest["outputs"].values():
        assert hashlib.sha256((actual_prefix.parent/item["file"]).read_bytes()).hexdigest()==item["sha256"]
    assert native_manifest["policy"].pop("TE_routing_policy") == "TE-routing-v1"
    assert native_manifest["policy"].pop("viral_evidence_scope").startswith("localized clusters,")
    assert native_manifest["policy"].pop("harmonizer_max_rank") == "species"
    assert native_manifest["policy"].pop("taxonomy_reporting_cap").startswith("reporting only,")
    assert native_manifest["policy"]["tesorter_retrovirus"].startswith("TE-routing-v1:")
    native_manifest["policy"]["tesorter_retrovirus"]=reference_manifest["policy"]["tesorter_retrovirus"]
    if rna and any(r["final_sequence_id"]=="retro" for r in baseline):
        reference_manifest["policy"]["rna_pair_homology_floor_demoted_ids"]=sorted(
            reference_manifest["policy"]["rna_pair_homology_floor_demoted_ids"]+["retro"])
    for label in ("taxonomy_reporting_audit","candidate_routing_audit","TE_candidates_fasta", "TE_candidates_summary",
                  "viral_candidates_fasta","viral_candidates_summary"):
        native_manifest["outputs"].pop(label)
    assert native_manifest.pop("reporting_bundle")["schema_version"] == "reporting-bundle-v1"
    for label in report_labels:
        native_manifest["outputs"].pop(label)
    reference_manifest["outputs"].pop("review_candidates_fasta")
    for manifest in (reference_manifest,native_manifest):
        for output in manifest["outputs"].values():
            output["file"] = Path(output["file"]).name
            output["sha256"]="contents_compared_above"
    assert native_manifest == reference_manifest


def run_both(args, reference, tmp_path, optimized_run=harmony.run):
    reference_dir = tmp_path / "reference"
    actual_dir = tmp_path / "indexed"
    reference_dir.mkdir()
    actual_dir.mkdir()
    baseline_args = Namespace(**vars(args))
    actual_args = Namespace(**vars(args))
    actual_args.max_rank = "species"  # Keep the frozen performance comparison at its original reporting depth.
    baseline_args.output_prefix = reference_dir / "sample"
    actual_args.output_prefix = actual_dir / "sample"
    reference.run(baseline_args)
    optimized_run(actual_args)
    compare_outputs(baseline_args.output_prefix, actual_args.output_prefix)


def run_cli(args, use_default_workers=False):
    command = [sys.executable, "-B", str(ROOT / "bin/run_viharmony.py")]
    for name, value in vars(args).items():
        if name == "workers" and use_default_workers:
            continue
        option = "--" + name.replace("_", "-")
        if isinstance(value, bool):
            if value:
                command.append(option)
        elif isinstance(value, list):
            command.extend([option, *map(str, value)])
        else:
            command.extend([option, str(value)])
    completed = subprocess.run(
        command, capture_output=True, text=True, timeout=60,
        env=dict(os.environ, PYTHONHASHSEED="17", PYTHONDONTWRITEBYTECODE="1"),
    )
    assert completed.returncode == 0, completed.stderr
    assert f"workers={1 if use_default_workers else args.workers}" in completed.stdout


def test_all_existing_policy_functions_are_unchanged(frozen_reference):
    _, source = frozen_reference
    def functions(text):
        return {node.name: ast.dump(node) for node in ast.parse(text).body
                if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name not in {"run", "parse_args"}}
    baseline = functions(source)
    current = functions((ROOT / "bin/run_viharmony.py").read_text(encoding="utf-8"))
    for name in baseline:
        assert current[name] == baseline[name], name


def test_extracted_candidate_policy_is_only_a_state_access_refactor(frozen_reference):
    """Check even unexercised branches of the old loop, not just its helpers."""
    _, source = frozen_reference
    original_run = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == "run")
    original_loop = next(n for n in original_run.body if isinstance(n, ast.For)
                         and ast.unparse(n.target) == "(final_id, sequence)")
    stop = next(i for i, n in enumerate(original_loop.body) if isinstance(n, ast.If)
                and ast.unparse(n.test) == "viral_decision in PRIMARY_VIRAL_DECISIONS")
    original = ast.Module(body=original_loop.body[:stop], type_ignores=[])
    current_function = next(n for n in ast.parse((ROOT / "bin/run_viharmony.py").read_text()).body
                            if isinstance(n, ast.FunctionDef) and n.name == "evaluate_candidate")
    start = next(i for i, n in enumerate(current_function.body) if isinstance(n, ast.Assign)
                 and ast.unparse(n.targets[0]) == "region")
    # This exact adapter is the approved TE-policy change, not a performance
    # refactor. Pin it and the reporting envelope; all other branches stay exact.
    reporting_start = next(i for i,n in enumerate(current_function.body)
        if isinstance(n,ast.Assign) and ast.unparse(n.targets[0]) == "reporting_audit")
    expected_reporting = ast.parse("reporting_audit = apply_taxonomy_reporting_cap(metadata, taxonomy, context.max_rank, context.sample_id)\nrouting_audit = te_routing.describe_candidate(metadata, evidence_audit_rows, reporting_audit, taxonomy_audit_rows, context.sample_id, previous_decision, previous_confidence)\nmetadata.update({name: routing_audit[name] for name in te_routing.METADATA_COLUMNS})\nif routing_audit['candidate_destination'] == 'TE_candidates':\n    metadata['viral_decision'] = 'TE_candidate'\nif review is not None:\n    review.update(metadata)").body
    assert ast.dump(ast.Module(body=current_function.body[reporting_start:-1],type_ignores=[])) == ast.dump(ast.Module(body=expected_reporting,type_ignores=[]))
    current = ast.Module(body=current_function.body[start:reporting_start],type_ignores=[])
    approved_origin = ast.parse("interpretation = tesorter_summary['sequence_interpretation']\nprevious_decision = adjudicate_viral_decision(discovery_status, record_type, strong_tools, qualified_tools, cellular, plasmid, interpretation)\nviral_decision = adjudicate_viral_decision(discovery_status, record_type, strong_tools, qualified_tools, cellular, plasmid, 'viral_candidate' if interpretation in te_routing.TE_VETO else interpretation)\nboundary_primary_eligible, provirus_retention_basis = adjudicate_disputed_ct3_boundary(region, rows, selected_vicat_rows)\nfloor_qualified_tools = rna_floor_qualified_tools(rows, qualified_tools)\nif pair_floor and context.input_type == 'rna' and (not context.has_locus_evidence) and (not strong_tools) and (floor_qualified_tools == {'deep6', 'vicat'}):\n    raise ValueError('RNA pair homology floor requires viCAT locus evidence')\n\ndef preserve_origin_safeguards(decision):\n    exception_loci = []\n    if not context.disable_rna_checkv_exception:\n        exception_loci = rna_checkv_exception_loci(context.input_type, decision, record_type, cellular, reported_plasmid_conflicts, strong_tools | qualified_tools, interpretation, task.loci, parent_id, task.length)\n    if exception_loci:\n        decision = 'retained_viral'\n    if decision == 'retained_provirus' and (not boundary_primary_eligible):\n        decision = 'provisional_provirus'\n    decision, floor_status = apply_rna_pair_homology_floor(decision, pair_floor, context.input_type, strong_tools, floor_qualified_tools, task.loci, parent_id, coords if record_type != 'input_contig' else '')\n    return (decision, floor_status, exception_loci)\nordinary_input = viral_decision\nviral_decision, pair_floor_status, checkv_exception_loci = preserve_origin_safeguards(ordinary_input)\nprevious_decision, previous_floor_status, _ = preserve_origin_safeguards(previous_decision) if previous_decision != ordinary_input else (viral_decision, pair_floor_status, checkv_exception_loci)\nprevious_confidence = 'provisional' if previous_floor_status == 'demoted' else viral_confidence\nif pair_floor_status == 'demoted':\n    viral_confidence = 'provisional'").body
    expressions = {
        "regions[final_id]": "task.region",
        "final_to_parent[final_id]": "task.parent_id",
        "headers[parent_id]": "task.header",
        "rows_by_final.get(final_id, [])": "task.evidence",
        "groups_by_sequence.get(final_id, [])": "task.groups",
        "len(sequence)": "task.length",
        "gate[parent_id]['discovery_status']": "task.discovery_status",
        "pair_loci": "task.loci",
        "args.sample_id": "context.sample_id",
        "args.input_type": "context.input_type",
        "getattr(args, 'disable_rna_checkv_exception', False)": "context.disable_rna_checkv_exception",
        "any(path.name.endswith('.vicat_orf_evidence.tsv') for path in args.evidence)": "context.has_locus_evidence",
        # The sole taxonomy call now passes through an optional taxonomy-only
        # adapter. Its original engine and all retention branches stay exact.
        "decide_taxonomy(rows, msl, minimum_vcontact3_length)": "decide_candidate_taxonomy(rows, msl, minimum_vcontact3_length, context, reconciliation_audit_rows)",
    }
    replacements = {ast.dump(ast.parse(old, mode="eval").body): ast.parse(new, mode="eval").body
                    for old, new in expressions.items()}

    class LiftState(ast.NodeTransformer):
        def visit(self, node):
            replacement = replacements.get(ast.dump(node))
            return copy.deepcopy(replacement) if replacement is not None else super().visit(node)

        def visit_Expr(self, node):
            call = node.value
            if (isinstance(call, ast.Call) and isinstance(call.func, ast.Attribute)
                    and call.func.attr == "append" and isinstance(call.func.value, ast.Name)):
                collection = call.func.value.id
                if collection in {"metadata_rows", "rna_checkv_exception_ids", "pair_floor_demoted_ids"}:
                    return None  # These collections are updated by the ordered parent.
                if collection == "review_rows":
                    return ast.Assign(targets=[ast.Name(id="review", ctx=ast.Store())], value=call.args[0])
            return self.generic_visit(node)

    transformed = LiftState().visit(original)
    current_a = next(i for i,n in enumerate(current.body) if isinstance(n,ast.Assign)
                     and ast.unparse(n.targets[0]) == "interpretation")
    current_b = next(i for i,n in enumerate(current.body) if isinstance(n,ast.Assign)
                     and ast.unparse(n.targets[0]) == "entity_interpretation")
    assert ast.dump(ast.Module(body=current.body[current_a:current_b],type_ignores=[])) == ast.dump(ast.Module(body=approved_origin,type_ignores=[]))
    original_a = next(i for i,n in enumerate(transformed.body) if isinstance(n,ast.Assign)
                      and ast.unparse(n.targets[0]) == "viral_decision")
    original_b = next(i for i,n in enumerate(transformed.body) if isinstance(n,ast.Assign)
                      and ast.unparse(n.targets[0]) == "entity_interpretation")
    current.body[current_a:current_b] = copy.deepcopy(transformed.body[original_a:original_b])
    assert ast.dump(transformed) == ast.dump(current)


@pytest.mark.parametrize("method", [
    "test_disputed_boundary_is_routed_to_provisional_outputs",
    "test_combines_parent_evidence_region_taxonomy_and_original_ids",
    "test_zero_candidates_still_writes_disposition_and_all_final_headers",
])
@pytest.mark.parametrize("workers", [1, 2, 4, 12])
def test_existing_end_to_end_cases_match_frozen_build(method, workers, frozen_reference, tmp_path):
    reference, _ = frozen_reference
    optimized_run = harmony.run
    calls = []
    def replay(args):
        args.workers = workers
        run_both(args, reference, tmp_path, optimized_run)
        calls.append(True)
        optimized_run(args)  # Preserve the fixture's own assertions as well.
    with patch.object(harmony, "run", side_effect=replay):
        getattr(existing_tests.ViharmonyTests(method), method)()
    assert calls == [True]


@pytest.mark.parametrize("input_type", ["dna", "rna"])
@pytest.mark.parametrize("audit_mode", ["none", "compact", "full"])
@pytest.mark.parametrize("workers,launch", [
    (1, "direct"), (2, "direct"), (4, "direct"), (12, "direct"),
    (1, "cli-default"), (12, "cli"),
])
def test_mixed_parent_loci_and_dispositions_match_frozen_build(input_type, audit_mode, workers, launch, frozen_reference, tmp_path):
    reference, _ = frozen_reference
    ids = ["exception", "floor_pass", "floor_fail", "nonviral", "ambiguous", "retro", *[f"extra_{i}" for i in range(48)], "unused"]
    normalized = tmp_path / "normalized.fasta"
    refined = tmp_path / "refined.fasta"
    normalized.write_text("".join(f">{key}\n{'A' * 2000}\n" for key in ids), encoding="utf-8")
    refined.write_text("".join(f">{key}\n{'A' * 2000}\n" for key in ids[:-1]), encoding="utf-8")
    header, gate, regions, evidence, loci, msl = [tmp_path / name for name in (
        "header.tsv", "gate.tsv", "regions.tsv", "evidence.tsv", "sample.vicat_orf_evidence.tsv", "MSL.csv")]
    write_tsv(header, ["sample_id", "sequence_id", "original_id", "original_header", "length"], [
        dict(sample_id="sample", sequence_id=key, original_id="shared-original", original_header=key, length="2000") for key in ids])
    write_tsv(gate, ["sequence_id", "advance_to_refinement", "discovery_status"], [
        dict(sequence_id=key, advance_to_refinement=str(key != "unused").lower(), discovery_status="ambiguous" if key in {"nonviral", "ambiguous"} else "viral") for key in ids])
    write_tsv(regions, ["sequence_id", "parent_sequence_id", "record_type", "coordinates", "original_length", "refined_length"], [
        dict(sequence_id=key, parent_sequence_id=key, record_type="input_contig", coordinates="", original_length="2000", refined_length="2000") for key in ids[:-1]])
    rows = []
    for key in ids[:-1]:
        tools = ["deep6", "virbot", "vicat"] if key == "exception" else ["deep6", "vicat"]
        if key == "nonviral":
            tools = ["vicat"]
        for tool in tools:
            rows.append(dict(sample_id="sample", sequence_id=key, parent_sequence_id="", tool=tool,
                             classification="virus", evidence_strength="qualified"))
        if key in {"exception", "nonviral", "ambiguous"}:
            rows.append(dict(sample_id="sample", sequence_id=key, parent_sequence_id="", tool="checkv",
                             classification="cellular", evidence_strength="qualified"))
        if key == "retro":
            rows.append(dict(sample_id="sample", sequence_id=key, parent_sequence_id="", tool="tesorter",
                             classification="retroelement", evidence_strength="strong", tesorter_evidence_category="retroelement"))
    write_tsv(evidence, ["sample_id", "sequence_id", "parent_sequence_id", "tool", "classification", "evidence_strength", "tesorter_evidence_category"], rows)
    locus_rows = [dict(sample_id="sample", sequence_id=key, locus_id=key + "_orf", coordinates="51-1950",
                      locus_classification="viral_supported", protein_length="630", best_bitscore="1100",
                      best_identity="95", best_query_coverage="95", best_subject_coverage="95")
                  for key in ("floor_pass", "exception", "unused", "floor_fail")]
    locus_rows[-1]["best_bitscore"] = "99"
    locus_rows.append(dict(locus_rows[0], sequence_id="unused", locus_classification="cellular_supported"))
    write_tsv(loci, list(locus_rows[0]), locus_rows)
    msl.write_text("Realm,Kingdom,Phylum,Class,Order,Family,Genus,Species\nRiboviria,,,,,,,\n", encoding="utf-8")
    args = Namespace(sample_id="sample", input_type=input_type, workers=workers, normalized_fasta=normalized,
                     header_map=header, discovery_gate=gate, refined_fasta=refined, region_map=regions,
                     ictv_msl=msl, evidence=[evidence, loci], vcontact3_groups=[], audit_mode=audit_mode,
                     rna_pair_homology_floor=True, output_prefix=tmp_path / "sample")
    optimized_run = harmony.run if launch == "direct" else (
        lambda value: run_cli(value, use_default_workers=launch == "cli-default")
    )
    run_both(args, reference, tmp_path, optimized_run)
    if input_type == "rna":
        with artifact(tmp_path / "indexed/sample.final_metadata.tsv").open(encoding="utf-8") as handle:
            metadata = {row["final_sequence_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        assert metadata["exception"]["rna_checkv_conflict_exception"] == "true"
        assert metadata["floor_pass"]["viral_decision"] == "retained_viral"
        assert metadata["floor_fail"]["viral_decision"] == "provisional_viral"
        assert metadata["nonviral"]["viral_decision"] == "likely_nonviral"
        assert metadata["ambiguous"]["viral_decision"] == "ambiguous_review"
