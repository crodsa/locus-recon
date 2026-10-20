#!/usr/bin/env python3
"""Compare fresh runs of the LIBA-6656 standard workflow with this deposit.

The graph and paths in this directory come from one assembly. This script
takes the output of new runs of the two commands in README.md, the standard
workflow and the terminal-path export, and records how they relate to the
deposit without editing it:

- whether each new local graph is the deposited graph: the same segment
  sequences, read on either strand, and the same links once segments are
  matched by sequence, so that segment numbering does not count; and which
  segments differ in their coverage tags;
- whether the exported terminal paths carry the same sequences under the
  same names;
- the tcdB sequence each run reported, with its tier, disposition and flags,
  and whether repeated runs agree.

Usage:
    python validation/liba6656-gfa-rerun/verify_release_rerun.py \\
        --run-dir run_4_threads/output --run-dir run_8_threads/output \\
        --terminal-paths tcdB_terminal_paths_8000_9000.fasta \\
        --output validation/liba6656-gfa-rerun/RELEASE_VERIFICATION.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

DEPOSIT = Path(__file__).resolve().parent
REPO_ROOT = DEPOSIT.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation.common import sha256_file, software_record, write_json  # noqa: E402

SAMPLE = "LIBA6656_ST154"
LOCUS = "tcdB"
GRAPH = "assembly_graph_after_simplification.gfa"
PATHS = "tcdB_terminal_paths_8000_9000.fasta"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FLIP = {"+": "-", "-": "+"}
PARAMETERS = ("threads", "memory_per_sample", "min_discovery_identity",
              "min_hsp_length", "region_padding", "cleanup")
REPORT_FIELDS = (
    "allele_length", "sequence_confidence", "result_disposition",
    "nearest_allele", "nearest_allele_identity_pct", "extracted_read_names",
    "remap_mean_depth", "mixed_site_count", "median_mixed_fraction",
    "max_alt_fraction", "strand_biased_sites", "span_clipped_bp",
    "span_clipped_start_bp", "span_clipped_end_bp", "span_continuation_contig",
    "span_continuation_bp", "qc_length_mod3",
)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def canonical_graph(path: Path) -> tuple[dict[str, list[str]], set[tuple]]:
    """Segment tags keyed by strand-independent sequence, and the links.

    A segment is named by the smaller of its sequence and that sequence's
    reverse complement, and a link by the smaller of its two equivalent
    spellings, so two assemblies of one graph compare equal whatever numbers
    the assembler gave their segments.
    """
    names: dict[str, tuple[str, bool]] = {}
    tags: dict[str, list[str]] = {}
    raw_links = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                sequence = fields[2].upper()
                key = min(sequence, reverse_complement(sequence))
                if key in tags:
                    raise ValueError(f"{path}: two segments share one sequence")
                names[fields[1]] = (key, key != sequence)
                tags[key] = sorted(fields[3:])
            elif fields[0] == "L":
                raw_links.append(fields[1:6])
    links = set()
    for source, source_strand, target, target_strand, overlap in raw_links:
        source_key, source_flipped = names[source]
        target_key, target_flipped = names[target]
        if source_flipped:
            source_strand = FLIP[source_strand]
        if target_flipped:
            target_strand = FLIP[target_strand]
        forward = (source_key, source_strand, target_key, target_strand, overlap)
        reverse = (target_key, FLIP[target_strand], source_key,
                   FLIP[source_strand], overlap)
        links.add(min(forward, reverse))
    return tags, links


def compare_graphs(deposited: Path, rerun: Path) -> dict:
    deposited_tags, deposited_links = canonical_graph(deposited)
    rerun_tags, rerun_links = canonical_graph(rerun)
    same_sequences = deposited_tags.keys() == rerun_tags.keys()
    shared = sorted(deposited_tags.keys() & rerun_tags.keys(), key=lambda key: (len(key), key))
    return {
        "segments": {"deposited": len(deposited_tags), "rerun": len(rerun_tags)},
        "links": {"deposited": len(deposited_links), "rerun": len(rerun_links)},
        "same_segment_sequences": same_sequences,
        "same_links_with_segments_matched_by_sequence": (
            same_sequences and deposited_links == rerun_links
        ),
        "coverage_tag_differences": [
            {"segment_length_bp": len(key), "deposited": deposited_tags[key],
             "rerun": rerun_tags[key]}
            for key in shared if deposited_tags[key] != rerun_tags[key]
        ],
    }


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name is not None:
                records[name].append(line.upper())
    return {key: "".join(value) for key, value in records.items()}


def compare_paths(deposited: Path, rerun: Path) -> dict:
    old, new = read_fasta(deposited), read_fasta(rerun)
    differing = sorted(name for name in old.keys() | new.keys()
                       if old.get(name) != new.get(name))
    return {
        "paths": {"deposited": len(old), "rerun": len(new)},
        "length_range_bp": {
            label: [min(map(len, records.values())), max(map(len, records.values()))]
            for label, records in (("deposited", old), ("rerun", new))
        },
        "same_sequence_under_same_name": not differing,
        "names_differing": differing,
    }


def split_flags(text: str) -> list[str]:
    """Split a qc_flags cell on the separators outside parentheses."""
    flags, depth, current = [], 0, []
    for character in text:
        depth += {"(": 1, ")": -1}.get(character, 0)
        if character == ";" and depth == 0:
            flags.append("".join(current).strip())
            current = []
        else:
            current.append(character)
    flags.append("".join(current).strip())
    return [flag for flag in flags if flag]


def deposited_checksums() -> dict[str, str]:
    with (DEPOSIT / "INPUT_CHECKSUMS.tsv").open(newline="") as handle:
        return {row["role"]: row["sha256"] for row in csv.DictReader(handle, delimiter="\t")}


def input_record(manifest: dict) -> dict:
    """SHA-256 of each run input beside the value this deposit records."""
    (sample,) = [entry for entry in manifest["samples"] if entry["sample_id"] == SAMPLE]
    recorded = deposited_checksums()
    record = {}
    for role, path in (("draft_assembly", sample["assembly"]["path"]),
                       ("read_1", sample["r1"]["path"]),
                       ("read_2", sample["r2"]["path"]),
                       ("ungapped_baits", manifest["bait"]["path"])):
        digest = sha256_file(path)
        record[role] = {"sha256": digest, "matches_deposit": digest == recorded[role]}
    return record


def run_record(run_dir: Path) -> tuple[dict, Path]:
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    with (run_dir / f"locus_recon_report_{LOCUS}.tsv").open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t")
                if row["sample_id"] == SAMPLE]
    if len(rows) != 1:
        raise ValueError(f"{run_dir}: expected one {SAMPLE} row, found {len(rows)}")
    sample_dir = run_dir / SAMPLE
    graph = sample_dir / "spades_local" / GRAPH
    reconstructed = sample_dir / f"{SAMPLE}_{LOCUS}_reconstructed.fasta"
    record = {
        "locus_recon_version": manifest["version"],
        "parameters": {key: manifest["parameters"][key] for key in PARAMETERS},
        "tools": {name: entry["version"] for name, entry in sorted(manifest["tools"].items())},
        "inputs": input_record(manifest),
        "graph_sha256": sha256_file(graph),
        "reconstructed_sha256": sha256_file(reconstructed),
        "report": {field: rows[0][field] for field in REPORT_FIELDS},
        "qc_flags": split_flags(rows[0]["qc_flags"]),
    }
    return record, graph


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--run-dir", type=Path, action="append", required=True,
                        help="main output directory of one standard-workflow run; "
                             "repeat for several runs")
    parser.add_argument("--terminal-paths", type=Path, required=True,
                        help="paths exported from the graph of the first run")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    records, graphs = zip(*(run_record(run_dir) for run_dir in args.run_dir))
    first = records[0]
    verification = {
        "schema_version": "1.0",
        "comparison": "fresh runs of the standard workflow and the terminal-path "
                      "export against this deposit",
        "software": software_record(REPO_ROOT),
        "runs": list(records),
        "runs_agree": all(
            record[key] == first[key]
            for record in records
            for key in ("graph_sha256", "reconstructed_sha256", "report", "qc_flags")
        ),
        "graph": compare_graphs(DEPOSIT / GRAPH, graphs[0]),
        "terminal_paths": compare_paths(DEPOSIT / PATHS, args.terminal_paths),
    }
    write_json(args.output, verification)
    graph = verification["graph"]
    inputs_match = all(entry["matches_deposit"]
                       for record in records for entry in record["inputs"].values())
    print(f"inputs match deposit {inputs_match}; "
          f"graph: same sequences {graph['same_segment_sequences']}, same links "
          f"{graph['same_links_with_segments_matched_by_sequence']}, "
          f"{len(graph['coverage_tag_differences'])} segment(s) with other coverage tags; "
          f"paths identical {verification['terminal_paths']['same_sequence_under_same_name']}; "
          f"runs agree {verification['runs_agree']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
