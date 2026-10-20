"""Graph evidence for scaffold placeholders inside a reconstruction."""

import random

from locus_recon.assembly_graph import (
    assess_placeholder_junctions,
    find_interior_n_runs,
    junction_support,
)


def _sequence(length: int, seed: int) -> str:
    return "".join(random.Random(seed).choices("ACGT", k=length))


def _gfa(tmp_path, segments, links, name="graph.gfa"):
    lines = ["H\tVN:Z:1.0"]
    for segment_id, sequence in segments:
        lines.append(f"S\t{segment_id}\t{sequence}\tDP:f:60.0")
    for left, right in links:
        lines.append(f"L\t{left}\t+\t{right}\t+\t0M")
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_find_interior_n_runs_ignores_terminal_runs():
    sequence = "NN" + "ACGT" * 10 + "N" * 5 + "ACGT" * 10 + "NNN"
    assert find_interior_n_runs(sequence, min_run_bp=1) == [(43, 47)]
    assert find_interior_n_runs(sequence, min_run_bp=6) == []


def test_junction_support_accepts_a_contiguous_spelling():
    left, right = _sequence(300, 1), _sequence(300, 2)
    support = junction_support([left + right], left[-40:], right[:40])
    assert support["assessed"] and support["supported"]
    assert support["paths_with_both_flanks"] == 1
    assert support["gaps_bp"] == [0]


def test_junction_support_rejects_an_alternative_spelling():
    left, right = _sequence(300, 3), _sequence(300, 4)
    filler = _sequence(60, 5)
    support = junction_support([left + filler + right], left[-40:], right[:40])
    assert support["assessed"] and not support["supported"]
    assert support["gaps_bp"] == [60]


def test_junction_support_finds_a_reverse_complemented_path():
    from locus_recon.io import reverse_complement

    left, right = _sequence(300, 6), _sequence(300, 7)
    support = junction_support(
        [reverse_complement(left + right)], left[-40:], right[:40]
    )
    assert support["supported"]


def test_placeholder_is_graph_supported(tmp_path):
    left, right = _sequence(300, 8), _sequence(300, 9)
    gfa = _gfa(tmp_path, [("1", left), ("2", right)], [("1", "2")])
    candidate = left + "N" * 100 + right
    result = assess_placeholder_junctions(candidate, gfa)
    assert result["verdict"] == "GRAPH_SUPPORTED"
    assert result["runs"] == 1 and result["total_bp"] == 100
    assert result["details"][0]["start"] == 301
    assert result["flags"][0].startswith("PLACEHOLDER_JUNCTION_GRAPH_SUPPORTED")


def test_placeholder_is_ambiguous_when_the_graph_spells_extra_sequence(tmp_path):
    left, right, middle = _sequence(300, 10), _sequence(300, 11), _sequence(150, 12)
    gfa = _gfa(
        tmp_path,
        [("1", left), ("2", middle), ("3", right)],
        [("1", "2"), ("2", "3")],
    )
    candidate = left + "N" * 100 + right
    result = assess_placeholder_junctions(candidate, gfa)
    assert result["verdict"] == "AMBIGUOUS"
    assert result["details"][0]["gaps_bp"] == [150]
    assert result["flags"][0].startswith("PLACEHOLDER_JUNCTION_AMBIGUOUS")


def test_candidate_without_placeholder_is_not_assessed(tmp_path):
    left = _sequence(300, 13)
    gfa = _gfa(tmp_path, [("1", left)], [])
    result = assess_placeholder_junctions(left, gfa)
    assert result["verdict"] == "NO_PLACEHOLDER"
    assert result["flags"] == []


def test_missing_graph_is_reported_rather_than_scored(tmp_path):
    left, right = _sequence(300, 14), _sequence(300, 15)
    candidate = left + "N" * 100 + right
    result = assess_placeholder_junctions(candidate, str(tmp_path / "absent.gfa"))
    assert result["verdict"] == "NOT_ASSESSED"
    assert result["flags"][0].startswith("PLACEHOLDER_JUNCTION_NOT_ASSESSED")
