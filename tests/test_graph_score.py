"""Competitive scoring of enumerated graph paths."""

from locus_recon.graph_score import (
    aggregate_depth,
    parse_idxstats,
    score_paths,
    write_scored_summary,
)

IDXSTATS = (
    "graph_path_1\t5000\t8000\t0\n"
    "graph_path_2\t5200\t4000\t0\n"
    "*\t0\t0\t120\n"
)


def _depth_text(entries):
    return "".join(
        f"{reference}\t{position}\t{depth}\n"
        for reference, position, depth in entries
    )


def test_parse_idxstats_drops_the_unmapped_row():
    stats = parse_idxstats(IDXSTATS)
    assert set(stats) == {"graph_path_1", "graph_path_2"}
    assert stats["graph_path_1"] == {"length": 5000, "mapped_reads": 8000}


def test_aggregate_depth_separates_covered_positions():
    text = _depth_text([("p1", 1, 30), ("p1", 2, 0), ("p1", 3, 20), ("p2", 1, 0)])
    aggregate = aggregate_depth(text)
    assert aggregate["p1"] == {
        "positions": 3, "covered_positions": 2, "total_depth": 50,
    }
    assert aggregate["p2"]["covered_positions"] == 0


def test_score_paths_ranks_by_uniquely_anchored_breadth():
    idxstats = {
        "covered": {"length": 4, "mapped_reads": 100},
        "patchy":  {"length": 4, "mapped_reads": 400},
    }
    depth = aggregate_depth(_depth_text(
        [("covered", i, 25) for i in range(1, 5)]
        + [("patchy", 1, 90), ("patchy", 2, 0), ("patchy", 3, 0), ("patchy", 4, 90)]
    ))
    rows = score_paths(idxstats, depth, graph_depths={"covered": 61.5})
    assert [row["candidate_id"] for row in rows] == ["covered", "patchy"]
    assert rows[0]["rank"] == 1
    assert rows[0]["unique_breadth_pct"] == 100.0
    assert rows[0]["unsupported_bp"] == 0
    assert rows[0]["mean_graph_depth"] == 61.5
    # the path with more mapped reads is not the better path if the reads are
    # piled on two positions
    assert rows[1]["mapped_reads"] > rows[0]["mapped_reads"]
    assert rows[1]["unsupported_bp"] == 2


def test_write_scored_summary_is_readable(tmp_path):
    rows = score_paths(parse_idxstats(IDXSTATS), {})
    for row in rows:
        row["nodes"] = "1+,2+"
    output = tmp_path / "paths.tsv"
    write_scored_summary(rows, str(output))
    lines = output.read_text().splitlines()
    assert lines[0].split("\t")[:3] == ["rank", "candidate_id", "length"]
    assert len(lines) == 3
    assert all(len(line.split("\t")) == 9 for line in lines)


def test_indistinguishable_paths_are_reported_as_unranked():
    """Every read fitting several paths means the order carries no preference."""
    from locus_recon.graph_score import ranking_is_informative

    idxstats = {
        "p1": {"length": 3, "mapped_reads": 1400},
        "p2": {"length": 3, "mapped_reads": 1480},
    }
    depth = aggregate_depth(_depth_text(
        [(ref, position, 0) for ref in ("p1", "p2") for position in (1, 2, 3)]
    ))
    rows = score_paths(idxstats, depth)
    assert all(row["unique_breadth_pct"] == 0.0 for row in rows)
    assert all(row["unsupported_bp"] == 3 for row in rows)
    assert not ranking_is_informative(rows)
    # and an informative case is distinguished from it
    better = score_paths(
        idxstats,
        aggregate_depth(_depth_text([("p1", position, 40) for position in (1, 2, 3)])),
    )
    assert ranking_is_informative(better)
