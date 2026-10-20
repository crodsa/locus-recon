from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from locus_recon.assembly_graph import parse_gfa
from locus_recon.copy_number import (
    GraphContextResult,
    GraphTargetHit,
    discover_graph_target_hits,
    infer_graph_contexts,
    reconcile_copy_number,
)
from locus_recon.depth_ratio import DepthRatioResult


def _write_gfa(tmp_path: Path, segments: dict[str, str], links: list[tuple]) -> Path:
    path = tmp_path / "graph.gfa"
    lines = ["H\tVN:Z:1.2"]
    lines.extend(f"S\t{name}\t{sequence}\tDP:f:20" for name, sequence in segments.items())
    lines.extend(
        f"L\t{source}\t{source_orientation}\t{target}\t{target_orientation}\t{overlap}M"
        for source, source_orientation, target, target_orientation, overlap in links
    )
    path.write_text("\n".join(lines) + "\n")
    return path


def _hit(
    segment_id: str = "T",
    query_start: int = 1,
    query_end: int = 100,
    subject_start: int = 1,
    subject_end: int = 100,
) -> GraphTargetHit:
    return GraphTargetHit(
        segment_id=segment_id,
        query_start=query_start,
        query_end=query_end,
        subject_start=subject_start,
        subject_end=subject_end,
        identity=100.0,
        bitscore=200.0,
    )


def _matched_context_graph(tmp_path: Path, count: int):
    segments = {"T": "ACGT" * 25}
    links = []
    left_sequences = ["AACCGG", "AAGGCC", "ACACGG"]
    right_sequences = ["CCGGAA", "GGAACC", "CAGTGC"]
    for index in range(count):
        left = f"L{index + 1}"
        right = f"R{index + 1}"
        segments[left] = left_sequences[index]
        segments[right] = right_sequences[index]
        links.append((left, "+", "T", "+", 0))
        links.append(("T", "+", right, "+", 0))
    return parse_gfa(str(_write_gfa(tmp_path, segments, links)))


@pytest.mark.parametrize("count", [1, 2, 3])
def test_infer_matching_graph_context_counts(tmp_path, count):
    graph = _matched_context_graph(tmp_path, count)

    result = infer_graph_contexts(
        graph,
        [_hit()],
        bait_length=100,
        min_context_bp=4,
        max_context_nodes=8,
        max_context_paths=64,
        max_context_bp=100,
    )

    assert result.status == "MATCHED_CONTEXT_COUNT"
    assert result.left_context_count == count
    assert result.right_context_count == count
    assert result.context_count == count
    assert result.context_lower_bound == count


def test_context_traversal_crosses_short_shared_node(tmp_path):
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {
            "T": "ACGT" * 25,
            "L1": "AACCGG",
            "L2": "AAGGCC",
            "S": "AT",
            "R1": "CCGGAA",
            "R2": "GGAACC",
        },
        [
            ("L1", "+", "T", "+", 0),
            ("L2", "+", "T", "+", 0),
            ("T", "+", "S", "+", 0),
            ("S", "+", "R1", "+", 0),
            ("S", "+", "R2", "+", 0),
        ],
    )))

    result = infer_graph_contexts(
        graph, [_hit()], bait_length=100, min_context_bp=6,
    )

    assert result.status == "MATCHED_CONTEXT_COUNT"
    assert (result.left_context_count, result.right_context_count) == (2, 2)


def test_default_context_span_uses_target_component_overlap_only(tmp_path):
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {
            "T": "A" * 200,
            "L": "C" * 200,
            "R": "G" * 200,
            "U": "T" * 1200,
            "V": "A" * 1200,
        },
        [
            ("L", "+", "T", "+", 60),
            ("T", "+", "R", "+", 60),
            ("U", "+", "V", "+", 1000),
        ],
    )))

    result = infer_graph_contexts(graph, [_hit(subject_start=51, subject_end=150)], bait_length=100)

    assert result.min_context_bp == 120
    assert result.status == "MATCHED_CONTEXT_COUNT"
    assert result.context_count == 1


def test_split_target_uses_query_outer_ends(tmp_path):
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {
            "T1": "ACGT" * 13,
            "T2": "TGCA" * 13,
            "L1": "AACCGG",
            "L2": "AAGGCC",
            "R1": "CCGGAA",
            "R2": "GGAACC",
        },
        [
            ("L1", "+", "T1", "+", 0),
            ("L2", "+", "T1", "+", 0),
            ("T1", "+", "T2", "+", 0),
            ("T2", "+", "R1", "+", 0),
            ("T2", "+", "R2", "+", 0),
        ],
    )))
    hits = [
        _hit("T1", 1, 50, 1, 50),
        _hit("T2", 51, 100, 1, 50),
    ]

    result = infer_graph_contexts(
        graph, hits, bait_length=100, min_context_bp=4,
    )

    assert result.status == "MATCHED_CONTEXT_COUNT"
    assert result.context_count == 2
    assert result.target_segments == ["T1", "T2"]


def test_distant_graph_branches_do_not_inflate_internal_single_copy_locus(tmp_path):
    target = "A" * 6000 + "ACGT" * 25 + "C" * 6000
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {
            "T": target,
            "L1": "AACCGG", "L2": "AAGGCC",
            "R1": "CCGGAA", "R2": "GGAACC",
        },
        [
            ("L1", "+", "T", "+", 0), ("L2", "+", "T", "+", 0),
            ("T", "+", "R1", "+", 0), ("T", "+", "R2", "+", 0),
        ],
    )))

    result = infer_graph_contexts(
        graph,
        [_hit(subject_start=6001, subject_end=6100)],
        bait_length=100,
        min_context_bp=4,
        max_context_bp=5000,
    )

    assert result.status == "MATCHED_CONTEXT_COUNT"
    assert result.context_count == 1
    assert "CONTEXT_RESOLVED_WITHIN_TARGET_SEGMENT" in result.flags


def test_reverse_complement_contexts_are_counted_once(tmp_path):
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {"T": "ACGT" * 25, "L": "AACCGG", "R1": "AAAACC", "R2": "GGTTTT"},
        [
            ("L", "+", "T", "+", 0),
            ("T", "+", "R1", "+", 0),
            ("T", "+", "R2", "+", 0),
        ],
    )))

    result = infer_graph_contexts(
        graph, [_hit()], bait_length=100, min_context_bp=6,
    )

    assert result.status == "MATCHED_CONTEXT_COUNT"
    assert result.context_count == 1


def test_one_sided_context_is_a_lower_bound(tmp_path):
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {"T": "ACGT" * 25, "L1": "AACCGG", "L2": "AAGGCC"},
        [("L1", "+", "T", "+", 0), ("L2", "+", "T", "+", 0)],
    )))

    result = infer_graph_contexts(
        graph, [_hit()], bait_length=100, min_context_bp=4,
    )

    assert result.status == "ONE_SIDED_LOWER_BOUND"
    assert result.context_count is None
    assert result.context_lower_bound == 2


def test_asymmetric_contexts_are_not_an_exact_call(tmp_path):
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {
            "T": "ACGT" * 25,
            "L1": "AACCGG", "L2": "AAGGCC",
            "R1": "CCGGAA", "R2": "GGAACC", "R3": "CAGTGC",
        },
        [
            ("L1", "+", "T", "+", 0), ("L2", "+", "T", "+", 0),
            ("T", "+", "R1", "+", 0), ("T", "+", "R2", "+", 0),
            ("T", "+", "R3", "+", 0),
        ],
    )))

    result = infer_graph_contexts(
        graph, [_hit()], bait_length=100, min_context_bp=4,
    )

    assert result.status == "ASYMMETRIC_CONTEXTS"
    assert result.context_count is None
    assert result.context_lower_bound == 2


def test_incomplete_bait_coverage_withholds_graph_count(tmp_path):
    graph = _matched_context_graph(tmp_path, 2)

    result = infer_graph_contexts(
        graph,
        [_hit(query_end=50, subject_end=50)],
        bait_length=100,
        min_query_coverage=0.90,
        min_context_bp=4,
    )

    assert result.status == "INSUFFICIENT_TARGET_COVERAGE"
    assert result.context_count is None
    assert result.bait_covered_bp == 50


def test_graph_discovery_requires_id_for_multirecord_bait(tmp_path):
    graph = _matched_context_graph(tmp_path, 1)
    bait = tmp_path / "baits.fasta"
    bait.write_text(">23S\n" + "A" * 100 + "\n>gyrB\n" + "C" * 100 + "\n")

    with pytest.raises(ValueError, match="--bait-id is required"):
        discover_graph_target_hits(graph, bait, min_hit_bp=20)


def test_graph_discovery_rejects_zero_minimum_hit_length(tmp_path):
    graph = _matched_context_graph(tmp_path, 1)
    bait = tmp_path / "bait.fasta"
    bait.write_text(">23S\n" + "A" * 100 + "\n")

    with pytest.raises(ValueError, match="min_hit_bp must be positive"):
        discover_graph_target_hits(graph, bait, min_hit_bp=0)


def test_graph_discovery_retains_deterministic_nonoverlapping_hits(tmp_path, monkeypatch):
    graph = _matched_context_graph(tmp_path, 1)
    bait = tmp_path / "baits.fasta"
    bait.write_text(">23S description\n" + "A" * 100 + "\n")

    def fake_run(command, **kwargs):
        assert Path(command[command.index("-query") + 1]).is_file()
        assert Path(command[command.index("-subject") + 1]).is_file()
        return SimpleNamespace(
            returncode=0,
            stdout=(
                "T\t99.0\t100\t1\t100\t1\t100\t200\n"
                "T\t98.0\t50\t1\t50\t1\t50\t80\n"
                "R1\t95.0\t10\t1\t10\t1\t10\t20\n"
            ),
            stderr="",
        )

    monkeypatch.setattr("locus_recon.copy_number.subprocess.run", fake_run)
    result = discover_graph_target_hits(
        graph, bait, bait_id="23S", min_identity=80, min_hit_bp=20,
    )

    assert result.bait_id == "23S"
    assert result.bait_length == 100
    assert len(result.hits) == 1
    assert result.hits[0].segment_id == "T"
    assert result.min_hit_bp == 20
    assert result.blastn_command[2] == "<selected-bait.fasta>"


def test_cycle_or_traversal_limit_is_indeterminate(tmp_path):
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {"T": "ACGT" * 25, "L": "AACCGG", "S1": "AT", "S2": "GC"},
        [
            ("L", "+", "T", "+", 0),
            ("T", "+", "S1", "+", 0),
            ("S1", "+", "S2", "+", 0),
            ("S2", "+", "S1", "+", 0),
        ],
    )))

    result = infer_graph_contexts(
        graph, [_hit()], bait_length=100, min_context_bp=20,
        max_context_nodes=3,
    )

    assert result.status == "TRAVERSAL_LIMIT_REACHED"
    assert result.context_count is None
    assert set(result.flags) & {"CYCLE_ENCOUNTERED", "MAX_CONTEXT_NODES_REACHED"}


def test_distal_cycle_downgrades_matching_contexts_to_a_lower_bound(tmp_path):
    # T branches right into R1, R2 and S1; the S1 branch enters a cycle and is
    # pruned before it resolves. The two sides therefore agree on 2 only
    # because the walk stopped, so 2 is a lower bound and not a count.
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {
            "T": "ACGT" * 25,
            "L1": "AACCGG", "L2": "AAGGCC",
            "R1": "CCGGAA", "R2": "GGAACC",
            "S1": "AT", "S2": "GC",
        },
        [
            ("L1", "+", "T", "+", 0), ("L2", "+", "T", "+", 0),
            ("T", "+", "R1", "+", 0), ("T", "+", "R2", "+", 0),
            ("T", "+", "S1", "+", 0),
            ("S1", "+", "S2", "+", 0), ("S2", "+", "S1", "+", 0),
        ],
    )))

    result = infer_graph_contexts(
        graph, [_hit()], bait_length=100, min_context_bp=6,
    )

    assert result.status == "TRAVERSAL_LIMIT_REACHED"
    assert result.context_count is None
    assert result.context_lower_bound == 2
    assert "CYCLE_ENCOUNTERED" in result.flags


def test_distal_node_limit_downgrades_matching_contexts_to_a_lower_bound(tmp_path):
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {
            "T": "ACGT" * 25,
            "L1": "AACCGG", "L2": "AAGGCC",
            "R1": "CCGGAA", "R2": "GGAACC",
            "S1": "AT", "S2": "GC", "S3": "AATT",
        },
        [
            ("L1", "+", "T", "+", 0), ("L2", "+", "T", "+", 0),
            ("T", "+", "R1", "+", 0), ("T", "+", "R2", "+", 0),
            ("T", "+", "S1", "+", 0),
            ("S1", "+", "S2", "+", 0), ("S2", "+", "S3", "+", 0),
        ],
    )))

    result = infer_graph_contexts(
        graph, [_hit()], bait_length=100, min_context_bp=6,
        max_context_nodes=1,
    )

    assert result.status == "TRAVERSAL_LIMIT_REACHED"
    assert result.context_count is None
    assert result.context_lower_bound == 2
    assert "MAX_CONTEXT_NODES_REACHED" in result.flags


def test_matching_contexts_stay_exact_when_the_traversal_is_exhaustive(tmp_path):
    # The guard on the two tests above must not fire on a clean locus: with no
    # pruning, matching counts remain an exact count.
    graph = parse_gfa(str(_write_gfa(
        tmp_path,
        {
            "T": "ACGT" * 25,
            "L1": "AACCGG", "L2": "AAGGCC",
            "R1": "CCGGAA", "R2": "GGAACC",
        },
        [
            ("L1", "+", "T", "+", 0), ("L2", "+", "T", "+", 0),
            ("T", "+", "R1", "+", 0), ("T", "+", "R2", "+", 0),
        ],
    )))

    result = infer_graph_contexts(
        graph, [_hit()], bait_length=100, min_context_bp=6,
    )

    assert result.flags == []
    assert result.status == "MATCHED_CONTEXT_COUNT"
    assert result.context_count == 2


def _depth(
    call: str,
    *,
    dosage: float | None,
    ratio: float,
    ci_low: float,
    ci_high: float,
    reliable: bool = True,
) -> DepthRatioResult:
    return DepthRatioResult(
        locus_regions=["contig:1-100"],
        locus_median_depth=ratio * 100,
        locus_median_depth_unfiltered=ratio * 100,
        backbone_median_depth=100.0,
        backbone_mad=5.0,
        ratio=ratio,
        ratio_ci_low=ci_low,
        ratio_ci_high=ci_high,
        assembly_copies=1,
        copies_estimate=dosage,
        ambiguity_index=1.0,
        gc_locus=None,
        single_copy_rejected=(call == "MULTICOPY_DEPTH"),
        reliable=reliable,
        assembly_region_count=1,
        dosage_estimate=dosage,
        dosage_status="ESTIMATED" if dosage is not None else "NOT_ESTIMATED",
        depth_call=call,
    )


def _graph_result(
    status: str,
    *,
    left: int,
    right: int,
    count: int | None,
    lower_bound: int | None,
) -> GraphContextResult:
    return GraphContextResult(
        target_segments=["T"],
        bait_length=100,
        bait_covered_bp=100,
        bait_coverage_fraction=1.0,
        left_context_count=left,
        right_context_count=right,
        context_count=count,
        context_lower_bound=lower_bound,
        status=status,
        flags=[],
        left_context_signatures=[],
        right_context_signatures=[],
        min_context_bp=100,
        max_context_nodes=8,
        max_context_paths=64,
        max_context_bp=5000,
    )


def test_two_graph_contexts_override_shifted_multicopy_depth_quantitation():
    depth = _depth(
        "MULTICOPY_DEPTH", dosage=3.221, ratio=3.221, ci_low=2.896, ci_high=3.546,
    )
    graph = _graph_result(
        "MATCHED_CONTEXT_COUNT", left=2, right=2, count=2, lower_bound=2,
    )

    result = reconcile_copy_number(depth, graph)

    assert result.copy_number_call == 2
    assert result.copy_number_kind == "INTEGER_CONTEXT_COUNT"
    # The depth interval (2.896, 3.546) excludes the graph count of 2, so the
    # two lines of evidence are not in consensus on the number. The graph count
    # still wins, under a method name that does not claim agreement.
    assert result.copy_number_method == "GRAPH_COUNT_OVER_DISCORDANT_DEPTH"
    assert result.copy_number_status == "SUPPORTED"
    assert "DEPTH_GRAPH_NUMERIC_DISCORDANCE" in result.consensus_flags
    assert depth.dosage_estimate == 3.221


def test_one_graph_context_preserves_tandem_depth_dosage_for_aphA1():
    depth = _depth(
        "MULTICOPY_DEPTH", dosage=4.2, ratio=4.2, ci_low=3.8, ci_high=4.6,
    )
    graph = _graph_result(
        "MATCHED_CONTEXT_COUNT", left=1, right=1, count=1, lower_bound=1,
    )

    result = reconcile_copy_number(depth, graph)

    assert result.copy_number_call == 4.2
    assert result.copy_number_kind == "MEAN_DEPTH_DOSAGE"
    assert result.copy_number_method == "DEPTH_TANDEM_COMPATIBLE"
    assert result.copy_number_status == "SUPPORTED"


def test_single_context_and_single_copy_compatible_depth_return_one():
    depth = _depth(
        "SINGLE_COPY_COMPATIBLE", dosage=0.98, ratio=0.98, ci_low=0.9, ci_high=1.1,
    )
    graph = _graph_result(
        "MATCHED_CONTEXT_COUNT", left=1, right=1, count=1, lower_bound=1,
    )

    result = reconcile_copy_number(depth, graph)

    assert result.copy_number_call == 1
    assert result.copy_number_kind == "INTEGER_CONTEXT_COUNT"
    assert result.copy_number_method == "GRAPH_DEPTH_CONSENSUS"


def test_graph_multicopy_and_single_copy_depth_are_a_conflict():
    depth = _depth(
        "SINGLE_COPY_COMPATIBLE", dosage=1.0, ratio=1.0, ci_low=0.9, ci_high=1.1,
    )
    graph = _graph_result(
        "MATCHED_CONTEXT_COUNT", left=2, right=2, count=2, lower_bound=2,
    )

    result = reconcile_copy_number(depth, graph)

    assert result.copy_number_call is None
    assert result.copy_number_status == "EVIDENCE_CONFLICT"
    assert "GRAPH_MULTICOPY_DEPTH_SINGLE_COPY_COMPATIBLE" in result.consensus_flags


def test_one_sided_graph_and_multicopy_depth_return_lower_bound():
    depth = _depth(
        "MULTICOPY_DEPTH", dosage=2.8, ratio=2.8, ci_low=2.5, ci_high=3.1,
    )
    graph = _graph_result(
        "ONE_SIDED_LOWER_BOUND", left=2, right=0, count=None, lower_bound=2,
    )

    result = reconcile_copy_number(depth, graph)

    assert result.copy_number_call is None
    assert result.copy_number_kind == "LOWER_BOUND"
    assert result.copy_number_method == "GRAPH_LOWER_BOUND"
    assert result.copy_number_status == "LOWER_BOUND"
    assert result.copy_number_lower_bound == 2


def test_missing_graph_uses_reliable_depth_only():
    depth = _depth(
        "MULTICOPY_DEPTH", dosage=2.4, ratio=2.4, ci_low=2.1, ci_high=2.7,
    )

    result = reconcile_copy_number(depth, None)

    assert result.copy_number_call == 2.4
    assert result.copy_number_kind == "MEAN_DEPTH_DOSAGE"
    assert result.copy_number_method == "DEPTH_ONLY"
    assert result.copy_number_status == "SUPPORTED"


def test_unreliable_depth_without_graph_is_indeterminate():
    depth = _depth(
        "DEPTH_INDETERMINATE", dosage=None, ratio=0.0, ci_low=0.0, ci_high=0.0,
        reliable=False,
    )

    result = reconcile_copy_number(depth, None)

    assert result.copy_number_call is None
    assert result.copy_number_kind == "NOT_ESTIMATED"
    assert result.copy_number_method == "NONE"
    assert result.copy_number_status == "INDETERMINATE"


def _low_uniqueness_multicopy(dosage: float = 4.2):
    depth = _depth(
        "MULTICOPY_DEPTH", dosage=dosage, ratio=dosage, ci_low=dosage - 0.4,
        ci_high=dosage + 0.4, reliable=False,
    )
    depth.dosage_status = "LOWER_BOUND"
    return depth


def test_low_mapping_uniqueness_multicopy_depth_is_a_floor_not_indeterminate():
    """The MAPQ filter only removes depth, so the multicopy verdict stands."""
    result = reconcile_copy_number(_low_uniqueness_multicopy(), None)

    assert result.copy_number_call is None
    assert result.copy_number_kind == "LOWER_BOUND"
    assert result.copy_number_method == "DEPTH_ONLY"
    assert result.copy_number_status == "LOWER_BOUND"
    assert result.copy_number_lower_bound == 4.2
    assert "DEPTH_LOWER_BOUND" in result.consensus_flags


def test_low_uniqueness_tandem_keeps_the_depth_floor_over_one_context():
    graph = _graph_result(
        "MATCHED_CONTEXT_COUNT", left=1, right=1, count=1, lower_bound=1,
    )
    result = reconcile_copy_number(_low_uniqueness_multicopy(4.2), graph)

    assert result.copy_number_kind == "LOWER_BOUND"
    assert result.copy_number_method == "DEPTH_ONLY"
    assert result.copy_number_lower_bound == 4.2
    assert "DEPTH_GRAPH_NUMERIC_DISCORDANCE" not in result.consensus_flags


def test_low_uniqueness_uses_a_larger_integer_graph_floor():
    graph = _graph_result(
        "MATCHED_CONTEXT_COUNT", left=3, right=3, count=3, lower_bound=3,
    )
    result = reconcile_copy_number(_low_uniqueness_multicopy(2.1), graph)

    assert result.copy_number_kind == "LOWER_BOUND"
    assert result.copy_number_method == "GRAPH_LOWER_BOUND"
    assert result.copy_number_lower_bound == 3
    assert "GRAPH_EVIDENCE_UNRESOLVED" not in result.consensus_flags


def test_low_uniqueness_single_copy_compatible_depth_stays_indeterminate():
    depth = _depth(
        "SINGLE_COPY_COMPATIBLE", dosage=0.6, ratio=0.6, ci_low=0.5, ci_high=0.7,
        reliable=False,
    )
    depth.dosage_status = "LOWER_BOUND"
    result = reconcile_copy_number(depth, None)

    assert result.copy_number_status == "INDETERMINATE"
    assert "DEPTH_EVIDENCE_UNRELIABLE" in result.consensus_flags
