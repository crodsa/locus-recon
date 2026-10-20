"""Public graph/depth copy-number consensus command."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess

from . import VERSION
from .assembly_graph import parse_gfa
from .copy_number import (
    GraphContextResult,
    discover_graph_target_hits,
    infer_graph_contexts,
    reconcile_copy_number,
)
from .depth_ratio import (
    CALL_THRESHOLD_DEFAULT,
    DepthRatioResult,
    estimate_locus_copy_number,
)


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("cannot be negative")
    return parsed


def _positive_float(value: str) -> float:
    parsed = float(value)
    if not parsed > 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def _fraction(value: str) -> float:
    parsed = float(value)
    if not 0 < parsed <= 1:
        raise argparse.ArgumentTypeError("must be in (0, 1]")
    return parsed


def _percentage(value: str) -> float:
    parsed = float(value)
    if not 0 < parsed <= 100:
        raise argparse.ArgumentTypeError("must be in (0, 100]")
    return parsed


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tool_version(executable: str) -> str:
    try:
        result = subprocess.run(
            [executable, "--version"], capture_output=True, text=True, timeout=15,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return f"unavailable: {exc}"
    text = (result.stdout.strip() or result.stderr.strip()).splitlines()
    return text[0] if text else f"unavailable: exit {result.returncode}"


def _parse_spans(parser: argparse.ArgumentParser, values: list[str], loci: list[str]):
    if not values:
        return None
    if len(values) != len(loci):
        parser.error("--gene-span must be given once per --locus, in the same order")
    spans = []
    for value in values:
        try:
            start, end = value.replace(",", "").split("-")
            spans.append((int(start), int(end)))
        except ValueError:
            parser.error(f"--gene-span expects START-END, got {value!r}")
    return spans


def _graph_fields(graph: GraphContextResult | None) -> dict:
    if graph is None:
        return {
            "graph_target_segments": [],
            "graph_bait_length": None,
            "graph_bait_covered_bp": None,
            "graph_bait_coverage_fraction": None,
            "left_context_count": None,
            "right_context_count": None,
            "graph_context_count": None,
            "graph_context_lower_bound": None,
            "graph_context_status": "NOT_EVALUATED",
            "graph_context_flags": [],
            "left_context_signatures": [],
            "right_context_signatures": [],
            "graph_min_context_bp": None,
            "graph_max_context_nodes": None,
            "graph_max_context_paths": None,
            "graph_max_context_bp": None,
        }
    return {
        "graph_target_segments": graph.target_segments,
        "graph_bait_length": graph.bait_length,
        "graph_bait_covered_bp": graph.bait_covered_bp,
        "graph_bait_coverage_fraction": graph.bait_coverage_fraction,
        "left_context_count": graph.left_context_count,
        "right_context_count": graph.right_context_count,
        "graph_context_count": graph.context_count,
        "graph_context_lower_bound": graph.context_lower_bound,
        "graph_context_status": graph.status,
        "graph_context_flags": graph.flags,
        "left_context_signatures": graph.left_context_signatures,
        "right_context_signatures": graph.right_context_signatures,
        "graph_min_context_bp": graph.min_context_bp,
        "graph_max_context_nodes": graph.max_context_nodes,
        "graph_max_context_paths": graph.max_context_paths,
        "graph_max_context_bp": graph.max_context_bp,
    }


def _combined_record(
    depth: DepthRatioResult,
    graph: GraphContextResult | None,
    consensus,
    *,
    provenance: dict,
) -> dict:
    record = {
        "schema_version": "1.0",
        "software_name": "locus-recon",
        "software_version": VERSION,
    }
    record.update(asdict(depth))
    record.update(_graph_fields(graph))
    record.update(consensus.as_dict())
    record["provenance"] = provenance
    return record


def _tsv_row(record: dict) -> dict:
    row = {key: value for key, value in record.items() if key != "provenance"}
    for key, value in list(row.items()):
        if isinstance(value, list):
            if value and isinstance(value[0], tuple):
                row[key] = ";".join(f"{start}-{end}" for start, end in value)
            else:
                row[key] = ";".join(str(item) for item in value)
    provenance = record["provenance"]
    row["input_sha256_json"] = json.dumps(
        provenance["input_sha256"], sort_keys=True, separators=(",", ":")
    )
    row["tool_versions_json"] = json.dumps(
        provenance["tool_versions"], sort_keys=True, separators=(",", ":")
    )
    row["parameters_json"] = json.dumps(
        provenance["parameters"], sort_keys=True, separators=(",", ":")
    )
    return row


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="locus-recon-copy-number",
        description=(
            "Combine unchanged depth multiplicity with optional SPAdes graph "
            "contexts without conflating continuous dosage and integer contexts."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--locus", required=True, action="append")
    parser.add_argument("--reference")
    parser.add_argument("--exclude", action="append", default=[])
    parser.add_argument("--gene-span", action="append", default=[])
    parser.add_argument("--gene-length", type=_positive_int)
    parser.add_argument("--min-bq", type=_nonnegative_int, default=20)
    parser.add_argument("--min-mq", type=_nonnegative_int, default=20)
    parser.add_argument(
        "--call-threshold", type=_positive_float, default=CALL_THRESHOLD_DEFAULT,
    )
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--n-boot", type=_positive_int, default=1000)
    parser.add_argument("--max-backbone-bootstrap-blocks", type=_positive_int, default=512)
    parser.add_argument("--samtools", default="samtools")

    parser.add_argument("--gfa")
    parser.add_argument("--bait")
    parser.add_argument("--bait-id")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--graph-min-identity", type=_percentage, default=80.0)
    parser.add_argument("--graph-min-hit-bp", type=_positive_int)
    parser.add_argument("--graph-min-query-coverage", type=_fraction, default=0.90)
    parser.add_argument("--graph-min-context-bp", type=_positive_int)
    parser.add_argument("--graph-max-context-nodes", type=_positive_int, default=8)
    parser.add_argument("--graph-max-context-paths", type=_positive_int, default=64)
    parser.add_argument("--graph-max-context-bp", type=_positive_int, default=5000)

    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-tsv")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(args.gfa) != bool(args.bait):
        parser.error("--gfa and --bait must be supplied together")
    if args.bait_id and not args.bait:
        parser.error("--bait-id requires --bait")
    spans = _parse_spans(parser, args.gene_span, args.locus)
    if spans and not args.gene_length:
        parser.error("--gene-span requires --gene-length")

    input_paths = {"bam": args.bam}
    for name in ("reference", "gfa", "bait"):
        value = getattr(args, name)
        if value:
            input_paths[name] = value
    for name, value in input_paths.items():
        if not Path(value).is_file():
            parser.error(f"{name} file was not found: {value}")
    output_json = Path(args.output_json)
    output_tsv = Path(args.output_tsv) if args.output_tsv else None
    if output_tsv and output_tsv.resolve() == output_json.resolve():
        parser.error("--output-json and --output-tsv must be different paths")

    depth_parameters = {
        "locus_regions": args.locus,
        "gene_spans": spans,
        "gene_length": args.gene_length,
        "reference_fasta": args.reference,
        "exclude_regions": args.exclude,
        "min_bq": args.min_bq,
        "min_mq": args.min_mq,
        "call_threshold": args.call_threshold,
        "seed": args.seed,
        "n_boot": args.n_boot,
        "max_backbone_bootstrap_blocks": args.max_backbone_bootstrap_blocks,
        "samtools": args.samtools,
    }
    depth = estimate_locus_copy_number(args.bam, **depth_parameters)

    graph_result = None
    discovery = None
    if args.gfa:
        graph = parse_gfa(args.gfa)
        discovery = discover_graph_target_hits(
            graph,
            args.bait,
            bait_id=args.bait_id,
            blastn=args.blastn,
            min_identity=args.graph_min_identity,
            min_hit_bp=args.graph_min_hit_bp,
        )
        graph_result = infer_graph_contexts(
            graph,
            discovery.hits,
            bait_length=discovery.bait_length,
            min_query_coverage=args.graph_min_query_coverage,
            min_context_bp=args.graph_min_context_bp,
            max_context_nodes=args.graph_max_context_nodes,
            max_context_paths=args.graph_max_context_paths,
            max_context_bp=args.graph_max_context_bp,
        )
    consensus = reconcile_copy_number(depth, graph_result)

    parameters = {
        "depth": {key: value for key, value in depth_parameters.items() if key != "samtools"},
        "graph": {
            "min_identity": args.graph_min_identity,
            "min_hit_bp": args.graph_min_hit_bp,
            "min_query_coverage": args.graph_min_query_coverage,
            "min_context_bp": args.graph_min_context_bp,
            "max_context_nodes": args.graph_max_context_nodes,
            "max_context_paths": args.graph_max_context_paths,
            "max_context_bp": args.graph_max_context_bp,
        } if args.gfa else None,
    }
    provenance = {
        "input_paths": {name: str(Path(value).resolve()) for name, value in input_paths.items()},
        "input_sha256": {name: _sha256(value) for name, value in input_paths.items()},
        "tool_versions": {
            "samtools": _tool_version(args.samtools),
            "blastn": _tool_version(args.blastn) if args.gfa else None,
        },
        "parameters": parameters,
        "graph_discovery": discovery.as_dict() if discovery else None,
    }
    record = _combined_record(depth, graph_result, consensus, provenance=provenance)

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    if output_tsv:
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        row = _tsv_row(record)
        with output_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(row), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerow(row)

    print(
        f"copy_number_call={consensus.copy_number_call} "
        f"kind={consensus.copy_number_kind} method={consensus.copy_number_method} "
        f"status={consensus.copy_number_status}"
    )
    print(f"JSON: {output_json}")
    if output_tsv:
        print(f"TSV: {output_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
