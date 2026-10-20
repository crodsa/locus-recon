import pytest

from locus_recon.assembly_graph import (
    enumerate_terminal_paths,
    parse_gfa,
    spell_path,
    write_paths_fasta,
)


def test_parse_and_enumerate_bubble_paths(tmp_path):
    graph_path = tmp_path / "bubble.gfa"
    graph_path.write_text(
        "\n".join(
            [
                "H\tVN:Z:1.2",
                "S\t1\tAAAC\tDP:f:10",
                "S\t2\tACGG\tDP:f:5",
                "S\t3\tACTT\tDP:f:6",
                "S\t4\tGGTA\tDP:f:10",
                "S\t5\tTTTA\tDP:f:10",
                "L\t1\t+\t2\t+\t2M",
                "L\t1\t+\t3\t+\t2M",
                "L\t2\t+\t4\t+\t2M",
                "L\t3\t+\t5\t+\t2M",
                "",
            ]
        )
    )

    graph = parse_gfa(str(graph_path))
    paths = enumerate_terminal_paths(graph)

    assert sorted(path.sequence for path in paths) == ["AAACGGTA", "AAACTTTA"]
    assert all(len(path.nodes) == 3 for path in paths)

    output = tmp_path / "paths.fasta"
    assert write_paths_fasta(paths, str(output), prefix="candidate") == 2
    assert output.read_text().count(">candidate_") == 2


def test_spell_path_rejects_nonmatching_sequence_overlap(tmp_path):
    graph_path = tmp_path / "bad.gfa"
    graph_path.write_text(
        "S\t1\tAAAA\n"
        "S\t2\tCCCC\n"
        "L\t1\t+\t2\t+\t2M\n"
    )
    graph = parse_gfa(str(graph_path))

    with pytest.raises(ValueError, match="overlap does not match"):
        spell_path(graph, (("1", "+"), ("2", "+")))


def test_parse_rejects_missing_embedded_sequence(tmp_path):
    graph_path = tmp_path / "missing.gfa"
    graph_path.write_text("S\t1\t*\n")

    with pytest.raises(ValueError, match="no embedded sequence"):
        parse_gfa(str(graph_path))


def test_zero_length_overlap_is_supported(tmp_path):
    graph_path = tmp_path / "zero.gfa"
    graph_path.write_text(
        "S\t1\tAAAA\n"
        "S\t2\tCCCC\n"
        "L\t1\t+\t2\t+\t0M\n"
    )

    paths = enumerate_terminal_paths(parse_gfa(str(graph_path)))

    assert len(paths) == 1
    assert paths[0].sequence == "AAAACCCC"


def test_parse_rejects_duplicate_segment_ids(tmp_path):
    graph_path = tmp_path / "duplicate.gfa"
    graph_path.write_text("S\t1\tAAAA\nS\t1\tCCCC\n")

    with pytest.raises(ValueError, match="Duplicate GFA segment ID"):
        parse_gfa(str(graph_path))
