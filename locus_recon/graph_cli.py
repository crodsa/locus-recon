"""CLI for exporting unresolved locus paths from a small SPAdes GFA graph."""

from __future__ import annotations

import argparse
import os
import re

from . import VERSION
from .assembly_graph import enumerate_terminal_paths, parse_gfa, write_paths_fasta


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("cannot be negative")
    return parsed


def _safe_prefix(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", value):
        raise argparse.ArgumentTypeError(
            "must start with a letter or number and contain only letters, "
            "numbers, dots, underscores, or hyphens"
        )
    return value


def write_path_summary(paths: list, output_path: str, prefix: str) -> None:
    """Write graph-path metadata without embedding full sequences."""

    with open(output_path, "w") as handle:
        handle.write("candidate_id\tlength\tmean_graph_depth\tnodes\n")
        for index, path in enumerate(paths, 1):
            nodes = ",".join(
                f"{segment_id}{orientation}"
                for segment_id, orientation in path.nodes
            )
            handle.write(
                f"{prefix}_{index}\t{len(path.sequence)}\t"
                f"{path.mean_depth:.3f}\t{nodes}\n"
            )


def main() -> None:
    """Enumerate bounded terminal paths from a local SPAdes assembly graph."""

    parser = argparse.ArgumentParser(
        prog="locus-recon-graph-paths",
        description=(
            "Export alternative terminal paths from a small SPAdes GFA graph. "
            "This research-mode command does not decide which path is a true "
            "allele; candidates require competitive read validation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument(
        "--gfa", required=True, metavar="GFA",
        help="SPAdes assembly_graph_after_simplification.gfa from a local assembly.",
    )
    parser.add_argument(
        "--output", required=True, metavar="FASTA",
        help="FASTA destination for the retained graph paths.",
    )
    parser.add_argument(
        "--summary", metavar="TSV",
        help="Path metadata TSV (default: OUTPUT.paths.tsv).",
    )
    parser.add_argument(
        "--min-length", type=_nonnegative_int, default=0, metavar="BP",
        help="Discard paths shorter than this length.",
    )
    parser.add_argument(
        "--max-length", type=_nonnegative_int, default=0, metavar="BP",
        help="Discard paths longer than this length; zero disables the limit.",
    )
    parser.add_argument(
        "--max-paths", type=_positive_int, default=10000, metavar="N",
        help="Safety limit for paths enumerated before length filtering.",
    )
    parser.add_argument(
        "--max-nodes", type=_positive_int, default=1000, metavar="N",
        help="Safety limit for segments in any one path.",
    )
    parser.add_argument(
        "--prefix", type=_safe_prefix, default="graph_path", metavar="TEXT",
        help="FASTA and TSV candidate identifier prefix.",
    )
    args = parser.parse_args()

    if args.max_length and args.max_length < args.min_length:
        parser.error("--max-length must be zero or at least --min-length")
    gfa_path = os.path.abspath(os.path.expanduser(args.gfa))
    if not os.path.isfile(gfa_path):
        parser.error(f"GFA file was not found: {gfa_path}")

    graph = parse_gfa(gfa_path)
    all_paths = enumerate_terminal_paths(
        graph,
        max_paths=args.max_paths,
        max_nodes=args.max_nodes,
    )
    paths = [
        path
        for path in all_paths
        if len(path.sequence) >= args.min_length
        and (not args.max_length or len(path.sequence) <= args.max_length)
    ]
    if not paths:
        parser.error(
            "No terminal graph paths passed the requested length filters."
        )

    output_path = os.path.abspath(os.path.expanduser(args.output))
    summary_path = os.path.abspath(
        os.path.expanduser(args.summary or output_path + ".paths.tsv")
    )
    for path in (output_path, summary_path):
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)

    write_paths_fasta(paths, output_path, prefix=args.prefix)
    write_path_summary(paths, summary_path, args.prefix)
    print(
        f"Exported {len(paths)} of {len(all_paths)} terminal graph paths.\n"
        f"FASTA: {output_path}\n"
        f"Summary: {summary_path}\n"
        "These are unresolved candidates. Validate copy number, path phase, "
        "and every reported base with competitive read mapping."
    )


if __name__ == "__main__":
    main()
