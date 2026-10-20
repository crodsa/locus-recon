from locus_recon.assembly_graph import GraphPath
from locus_recon.graph_cli import write_path_summary


def test_write_path_summary(tmp_path):
    output = tmp_path / "paths.tsv"
    paths = [
        GraphPath(
            nodes=(("1", "+"), ("2", "-")),
            sequence="AACCGG",
            mean_depth=12.3456,
        )
    ]

    write_path_summary(paths, str(output), "candidate")

    assert output.read_text().splitlines() == [
        "candidate_id\tlength\tmean_graph_depth\tnodes",
        "candidate_1\t6\t12.346\t1+,2-",
    ]
