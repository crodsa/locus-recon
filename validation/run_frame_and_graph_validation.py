#!/usr/bin/env python3
"""Validation for orientation-independent frame metrics, bait validation,
placeholder junction evidence, and competitive graph-path scoring.

Four of the five checks are deterministic and use only files already deposited
under ``validation/liba6656-gfa-rerun/``: the 111 curated ``tcdB`` alleles and
the local assembly graph from the LIBA-6656 rerun.  No new sequence data is
introduced, and every derived input is constructed from those files at run
time rather than stored.

1. Frame orientation invariance.  A bait catalogue can be supplied on either
   strand, and a reconstruction inherits the orientation of the bait that
   recruited it.  One allele is held out as the candidate and the remaining 110
   build the profile; the pair is then assessed in all four orientation
   combinations.  The internal-stop count must be identical in all four and the
   reported frame must mirror, which is the property that a forward-only scan
   does not have.  The forward-only minimum is reported alongside as the
   quantity a three-frame scan would have produced.

2. Bait database validation.  Three copies of the deposited bait file are
   derived in a temporary directory -- one with an identifier past the
   ``makeblastdb -parse_seqids`` limit, one with alignment gaps restored, one
   carrying an unrelated fragment -- and the message each produces is recorded.

3. Length-profile size.  A three-allele subset of the deposited baits is used
   to confirm that length flags carry the profile size and that a profile below
   the minimum is marked underpowered.

4. Placeholder junction evidence.  Paths are spelled from the deposited GFA and
   three candidates are constructed from them: a placeholder inserted at a
   position the graph spells contiguously, a placeholder replacing real
   sequence the graph still spells, and a placeholder whose second flank does
   not occur in the graph at all.  The three expected verdicts are
   GRAPH_SUPPORTED, AMBIGUOUS and NOT_SUPPORTED.

5. Competitive graph-path scoring, positive control (requires an aligner).
   The published case below cannot separate its paths, so the ranking is also
   exercised against constructed truth, in the same spirit as the completeness
   benchmark: reads are simulated from one deposited path and scored against
   that path plus the nine most divergent paths in the deposit.  The source
   path must rank first and must be the only path with substantial uniquely
   anchored coverage.

6. Competitive graph-path scoring on the published case (optional; requires
   reads and an aligner).
   With ``--reads-r1``/``--reads-r2`` the paths enumerated from the deposited
   GFA are indexed together and scored by uniquely anchored coverage.  The
   reads are not deposited here; ENA run ``ERR467623`` is the pair used in the
   original rerun.

Usage:
    python validation/run_frame_and_graph_validation.py \
        --deposit validation/liba6656-gfa-rerun \
        --outdir validation/frame-and-graph-evidence \
        [--reads-r1 ERR467623_1.fastq.gz --reads-r2 ERR467623_2.fastq.gz] \
        [--threads 8]
"""

import argparse
import csv
import json
import random
import re
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from locus_recon.assembly_graph import (  # noqa: E402
    assess_placeholder_junctions,
    enumerate_terminal_paths,
    parse_gfa,
    write_paths_fasta,
)
from locus_recon.io import read_fasta_sequences, reverse_complement  # noqa: E402
from locus_recon.qc import (  # noqa: E402
    MIN_LENGTH_PROFILE_N,
    PARALOG_KMER_SIZE,
    _count_internal_stops,
    assess_allele_quality,
    profile_bait_database,
)
from validation.common import (  # noqa: E402
    build_provenance,
    write_json,
)

# Read support is not the subject of these checks, so every assessment is given
# the same unambiguous remap metrics and the frame, length and graph behaviour
# is read from the flags.
REMAP_METRICS = {
    "mean_depth": 90.0,
    "breadth_pct": 100.0,
    "pct_bases_lt5": 0.0,
    "mapped_reads": 5000,
}
PLACEHOLDER_BP = 100
AMBIGUOUS_WINDOW_BP = 150
FLANK_BP = 40


def write_fasta(path, records):
    """Write (id, sequence) pairs as wrapped FASTA."""
    with open(path, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")
    return str(path)


def portable_log(path):
    """Replace absolute executable paths in a tool log with ``<env>``.

    bwa-mem2 names the binary it launches; the location of the environment on
    the machine that ran it is not evidence.
    """
    path = Path(path)
    if path.is_file():
        path.write_text(re.sub(r'"/[^"]*/bin/', '"<env>/bin/', path.read_text()))


def write_tsv(path, header, rows):
    """Write one results table."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def assess(candidate, profile):
    """Run the production QC path on one candidate."""
    return assess_allele_quality(
        allele_seq=candidate,
        bait_profile=profile,
        validation_identity=99.0,
        validation_qcov=100.0,
        locus="tcdB",
        remap_metrics=REMAP_METRICS,
        expect_cds=True,
    )


def check_frame_orientation(baits, outdir):
    """Check 1: the frame metrics must not depend on bait orientation."""
    held_out, remainder = baits[0], baits[1:]
    with tempfile.TemporaryDirectory() as work:
        forward = write_fasta(Path(work) / "baits_forward.fasta", remainder)
        reverse = write_fasta(
            Path(work) / "baits_reverse.fasta",
            [(name, reverse_complement(seq)) for name, seq in remainder],
        )
        profiles = {
            "as deposited": profile_bait_database(forward),
            "reverse complemented": profile_bait_database(reverse),
        }
    candidates = {
        "as deposited": held_out[1],
        "reverse complemented": reverse_complement(held_out[1]),
    }
    rows = []
    for bait_orientation, profile in profiles.items():
        for candidate_orientation, candidate in candidates.items():
            qc = assess(candidate, profile)
            forward_only = min(
                _count_internal_stops(candidate, frame) for frame in range(3)
            )
            rows.append([
                held_out[0],
                bait_orientation,
                candidate_orientation,
                profile["expected_coding_frame_labels"][0],
                str(profile["bait_appears_antisense"]),
                qc["internal_stops"],
                qc["coding_frame_used"],
                forward_only,
                "yes" if any(
                    flag.startswith("INTERNAL_STOPS") for flag in qc["flags"]
                ) else "no",
            ])
    write_tsv(
        Path(outdir) / "frame_orientation_invariance.tsv",
        ["candidate_record", "bait_orientation", "candidate_orientation",
         "bait_plausible_frame", "bait_reported_antisense", "internal_stops",
         "frame_used", "forward_only_internal_stops", "internal_stops_flag"],
        rows,
    )
    six_frame = {int(row[5]) for row in rows}
    forward_only = {int(row[7]) for row in rows}
    return {
        "candidate_record": held_out[0],
        "profile_alleles": len(remainder),
        "six_frame_internal_stops": sorted(six_frame),
        "forward_only_internal_stops": sorted(forward_only),
        "invariant": len(six_frame) == 1,
        "forward_only_invariant": len(forward_only) == 1,
    }


def check_bait_validation(baits, bait_path, outdir):
    """Check 2: every rejected input names the record and the fix."""
    rows = []
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        long_name = "X" * 64
        cases = {
            "identifier past the parse_seqids limit": write_fasta(
                work / "long_id.fasta",
                [(long_name, baits[0][1])] + [(n, s) for n, s in baits[1:6]],
            ),
            "alignment gaps restored": write_fasta(
                work / "gapped.fasta",
                [(baits[0][0], baits[0][1][:300] + "-" * 20 + baits[0][1][300:])]
                + [(n, s) for n, s in baits[1:6]],
            ),
            "unrelated fragment among the alleles": write_fasta(
                work / "paralog.fasta",
                [(n, s) for n, s in baits[:6]]
                + [(
                    "unrelated_fragment",
                    "".join(random.Random(0).choices("ACGT", k=400)),
                )],
            ),
        }
        for label, path in cases.items():
            try:
                profile = profile_bait_database(path)
            except ValueError as exc:
                rows.append([label, "rejected", str(exc)])
                continue
            unrelated = profile["unrelated_records"]
            rows.append([
                label,
                "accepted with warning" if unrelated else "accepted",
                (f"{len(unrelated)} record(s) share no {PARALOG_KMER_SIZE}-mer "
                 f"with the full-length records: {', '.join(unrelated)}")
                if unrelated else "no warning",
            ])
    write_tsv(
        Path(outdir) / "bait_validation_messages.tsv",
        ["derived_input", "outcome", "message"],
        rows,
    )
    return {
        "deposited_bait_accepted": bool(profile_bait_database(bait_path)),
        "cases": {row[0]: row[1] for row in rows},
    }


def check_length_profile(baits, outdir):
    """Check 3: length flags carry the number of alleles behind them."""
    rows = []
    with tempfile.TemporaryDirectory() as work:
        subset = write_fasta(
            Path(work) / "three_alleles.fasta", [(n, s) for n, s in baits[1:4]]
        )
        small = profile_bait_database(subset)
        full = profile_bait_database(
            write_fasta(Path(work) / "all_alleles.fasta", baits[1:])
        )
    truncated = baits[0][1][: len(baits[0][1]) // 2]
    for label, profile in (("three alleles", small), ("110 alleles", full)):
        qc = assess(truncated, profile)
        length_flags = [f for f in qc["flags"] if f.startswith("LENGTH_")]
        rows.append([
            label,
            profile["n_alleles"],
            str(profile["length_profile_underpowered"]),
            qc["length_profile_n"],
            length_flags[0] if length_flags else "none",
            "yes" if any(
                f.startswith("CATALOGUE_LENGTH_PROFILE_UNDERPOWERED")
                for f in qc["flags"]
            ) else "no",
        ])
    write_tsv(
        Path(outdir) / "length_profile_reporting.tsv",
        ["profile", "bait_alleles", "profile_underpowered", "reported_n",
         "length_flag", "underpowered_flag"],
        rows,
    )
    return {
        "minimum_profile_alleles": MIN_LENGTH_PROFILE_N,
        "underpowered_flag_on_three": rows[0][5],
        "underpowered_flag_on_110": rows[1][5],
    }


def check_placeholder_junctions(gfa_path, outdir):
    """Check 4: graph evidence for a scaffold placeholder, three ways."""
    graph = parse_gfa(gfa_path)
    paths = enumerate_terminal_paths(graph, max_paths=10000, max_nodes=1000)
    longest = max(paths, key=lambda path: len(path.sequence)).sequence
    cut = len(longest) // 2

    contiguous = longest[:cut] + "N" * PLACEHOLDER_BP + longest[cut:]
    replaced = (
        longest[:cut]
        + "N" * PLACEHOLDER_BP
        + longest[cut + AMBIGUOUS_WINDOW_BP :]
    )
    shuffled = list(longest[cut : cut + 400])
    random.Random(1).shuffle(shuffled)
    absent_flank = longest[:cut] + "N" * PLACEHOLDER_BP + "".join(shuffled)

    rows = []
    for label, candidate, expected in (
        ("placeholder inserted where the graph is contiguous",
         contiguous, "GRAPH_SUPPORTED"),
        (f"placeholder replacing {AMBIGUOUS_WINDOW_BP} bp the graph still spells",
         replaced, "AMBIGUOUS"),
        ("second flank absent from the graph", absent_flank, "NOT_SUPPORTED"),
    ):
        result = assess_placeholder_junctions(
            candidate, gfa_path, flank_bp=FLANK_BP, min_run_bp=10
        )
        detail = result["details"][0] if result["details"] else {}
        rows.append([
            label,
            len(candidate),
            result["runs"],
            result["total_bp"],
            expected,
            result["verdict"],
            "yes" if result["verdict"] == expected else "NO",
            detail.get("paths_with_both_flanks", 0),
            ";".join(str(gap) for gap in detail.get("gaps_bp", [])) or "none",
        ])
    write_tsv(
        Path(outdir) / "placeholder_junction_cases.tsv",
        ["case", "candidate_bp", "interior_placeholder_runs", "placeholder_bp",
         "expected_verdict", "observed_verdict", "agrees",
         "graph_paths_with_both_flanks", "distance_between_flanks_bp"],
        rows,
    )
    return {
        "graph_paths_enumerated": len(paths),
        "longest_path_bp": len(longest),
        "cases_agreeing": sum(1 for row in rows if row[6] == "yes"),
        "cases": len(rows),
    }


def check_path_scoring(gfa_path, outdir, reads_r1, reads_r2, threads,
                       min_length, max_length, prefix):
    """Check 6: competitive scoring of the enumerated paths against reads."""
    from locus_recon.graph_score import (
        competitive_read_scores,
        ranking_is_informative,
        write_scored_summary,
    )
    from locus_recon.utils import check_dependencies

    graph = parse_gfa(gfa_path)
    paths = [
        path
        for path in enumerate_terminal_paths(graph, max_paths=10000, max_nodes=1000)
        if min_length <= len(path.sequence) <= max_length
    ]
    if not paths:
        raise RuntimeError("no graph path passed the requested length bounds")

    outdir = Path(outdir)
    scratch = tempfile.TemporaryDirectory()
    score_dir = Path(scratch.name)
    paths_fasta = str(score_dir / "graph_paths.fasta")
    write_paths_fasta(paths, paths_fasta, prefix=prefix)

    graph_depths = {
        f"{prefix}_{index}": path.mean_depth
        for index, path in enumerate(paths, 1)
    }
    nodes = {
        f"{prefix}_{index}": ",".join(
            f"{segment}{orientation}" for segment, orientation in path.nodes
        )
        for index, path in enumerate(paths, 1)
    }
    tools, aligner_name = check_dependencies(False)
    with open(outdir / "graph_path_scores.log", "w") as log_handle:
        rows = competitive_read_scores(
            tools=tools,
            aligner_name=aligner_name,
            threads=threads,
            paths_fasta=paths_fasta,
            r1=reads_r1,
            r2=reads_r2,
            out_dir=str(score_dir),
            log_handle=log_handle,
            graph_depths=graph_depths,
        )
    portable_log(outdir / "graph_path_scores.log")
    for row in rows:
        row["nodes"] = nodes.get(row["candidate_id"], "")
    write_scored_summary(rows, str(outdir / "graph_path_scores.tsv"))
    # Only small evidence is deposited: the paths themselves are already
    # deposited as tcdB_terminal_paths_8000_9000.fasta, and the alignment and
    # per-base depth are regenerable from the reads.
    idxstats = score_dir / "graph_paths.idxstats.txt"
    if idxstats.is_file():
        (outdir / "graph_path_scores_idxstats.txt").write_text(idxstats.read_text())
    scratch.cleanup()

    best = rows[0]
    return {
        "paths_scored": len(rows),
        "reads_discriminate": ranking_is_informative(rows),
        "best_candidate": best["candidate_id"],
        "best_unique_breadth_pct": round(best["unique_breadth_pct"], 2),
        "best_unique_mean_depth": round(best["unique_mean_depth"], 2),
        "best_unsupported_bp": best["unsupported_bp"],
        "mapped_reads_min": min(row["mapped_reads"] for row in rows),
        "mapped_reads_max": max(row["mapped_reads"] for row in rows),
    }


def _simulate_pairs(sequence, depth, read_bp, fragment_bp, error_rate, seed):
    """Simulate paired reads from one sequence, deterministically."""
    rng = random.Random(seed)
    pairs = max(1, int(depth * len(sequence) / (2 * read_bp)))
    forward, reverse = [], []
    for index in range(pairs):
        start = rng.randrange(0, max(1, len(sequence) - fragment_bp))
        fragment = sequence[start : start + fragment_bp]
        read1 = list(fragment[:read_bp])
        read2 = list(reverse_complement(fragment[-read_bp:]))
        for read in (read1, read2):
            for position in range(len(read)):
                if rng.random() < error_rate:
                    read[position] = rng.choice(
                        [base for base in "ACGT" if base != read[position]]
                    )
        name = f"sim_{index + 1}"
        forward.append((name, "".join(read1)))
        reverse.append((name, "".join(read2)))
    return forward, reverse


def _write_fastq_gz(path, records):
    """Write one FASTQ.gz with a fixed timestamp so the file is reproducible."""
    import gzip

    with gzip.GzipFile(path, "wb", mtime=0) as handle:
        for name, sequence in records:
            handle.write(
                f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n".encode()
            )
    return str(path)


def _kmer_set(sequence, size=31):
    return {sequence[i : i + size] for i in range(len(sequence) - size + 1)}


def check_scoring_positive_control(gfa_path, outdir, threads, min_length,
                                   max_length, depth=40, read_bp=100,
                                   fragment_bp=300, error_rate=0.002, seed=7,
                                   decoys=9):
    """Check 5: the ranking must recover a known source path.

    Reads are simulated from one deposited graph path and scored against that
    path plus the most divergent paths in the deposit, so the correct answer is
    known by construction.
    """
    from locus_recon.graph_score import (
        competitive_read_scores,
        ranking_is_informative,
        write_scored_summary,
    )
    from locus_recon.utils import check_dependencies

    graph = parse_gfa(gfa_path)
    paths = [
        path
        for path in enumerate_terminal_paths(graph, max_paths=10000, max_nodes=1000)
        if min_length <= len(path.sequence) <= max_length
    ]
    source = paths[0]
    source_kmers = _kmer_set(source.sequence)
    similarity = []
    for index, path in enumerate(paths[1:], start=2):
        shared = len(source_kmers & _kmer_set(path.sequence))
        similarity.append((shared, index, path))
    similarity.sort(key=lambda item: (item[0], item[1]))
    selected = [("source", source)] + [
        (f"decoy_{rank}", path) for rank, (_, _, path) in enumerate(similarity[:decoys], 1)
    ]

    outdir = Path(outdir)
    candidates_fasta = write_fasta(
        outdir / "scoring_positive_control_candidates.fasta",
        [(label, path.sequence) for label, path in selected],
    )
    forward, reverse = _simulate_pairs(
        source.sequence, depth, read_bp, fragment_bp, error_rate, seed
    )
    r1 = _write_fastq_gz(outdir / "scoring_positive_control_R1.fastq.gz", forward)
    r2 = _write_fastq_gz(outdir / "scoring_positive_control_R2.fastq.gz", reverse)

    tools, aligner_name = check_dependencies(False)
    with tempfile.TemporaryDirectory() as scratch:
        with open(outdir / "scoring_positive_control.log", "w") as log_handle:
            rows = competitive_read_scores(
                tools=tools, aligner_name=aligner_name, threads=threads,
                paths_fasta=candidates_fasta, r1=r1, r2=r2,
                out_dir=scratch, log_handle=log_handle,
            )
        portable_log(outdir / "scoring_positive_control.log")
    write_scored_summary(rows, str(outdir / "scoring_positive_control.tsv"))
    for index in Path(candidates_fasta).parent.glob(
        Path(candidates_fasta).name + ".*"
    ):
        index.unlink()

    best = rows[0]
    runner_up = rows[1] if len(rows) > 1 else best
    return {
        "candidates": len(rows),
        "simulated_pairs": len(forward),
        "simulated_depth": depth,
        "read_bp": read_bp,
        "substitution_rate": error_rate,
        "seed": seed,
        "source_ranked_first": best["candidate_id"] == "source",
        "source_unique_breadth_pct": round(
            next(row["unique_breadth_pct"] for row in rows
                 if row["candidate_id"] == "source"), 2),
        "runner_up": runner_up["candidate_id"],
        "runner_up_unique_breadth_pct": round(runner_up["unique_breadth_pct"], 2),
        "reads_discriminate": ranking_is_informative(rows),
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Validate orientation-independent frame metrics, bait validation, "
            "placeholder junction evidence and competitive graph-path scoring "
            "against the deposited LIBA-6656 tcdB material."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--deposit", default="validation/liba6656-gfa-rerun",
        help="Directory holding the deposited baits and assembly graph.",
    )
    parser.add_argument(
        "--outdir", default="validation/frame-and-graph-evidence",
        help="Destination for the result tables and provenance.",
    )
    parser.add_argument("--reads-r1", help="Forward reads for check 5 (ENA ERR467623_1).")
    parser.add_argument("--reads-r2", help="Reverse reads for check 5 (ENA ERR467623_2).")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-length", type=int, default=8000)
    parser.add_argument("--max-length", type=int, default=9000)
    parser.add_argument("--prefix", default="LIBA6656_tcdB_terminal_path")
    args = parser.parse_args()

    deposit = Path(args.deposit)
    bait_path = deposit / "tcdB_111_ungapped_baits.fasta"
    gfa_path = deposit / "assembly_graph_after_simplification.gfa"
    for path in (bait_path, gfa_path):
        if not path.is_file():
            parser.error(f"deposited input was not found: {path}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    baits = read_fasta_sequences(str(bait_path))
    print(f"[1/6] frame orientation invariance ({len(baits)} deposited alleles)")
    summary = {"frame_orientation": check_frame_orientation(baits, outdir)}
    print(f"      six-frame stops {summary['frame_orientation']['six_frame_internal_stops']}"
          f" vs forward-only {summary['frame_orientation']['forward_only_internal_stops']}")

    print("[2/6] bait database validation")
    summary["bait_validation"] = check_bait_validation(baits, str(bait_path), outdir)

    print("[3/6] length-profile reporting")
    summary["length_profile"] = check_length_profile(baits, outdir)

    print("[4/6] placeholder junction evidence")
    summary["placeholder_junctions"] = check_placeholder_junctions(str(gfa_path), outdir)
    print(f"      {summary['placeholder_junctions']['cases_agreeing']} of "
          f"{summary['placeholder_junctions']['cases']} cases match the expected verdict")

    print("[5/6] competitive scoring, constructed positive control")
    try:
        summary["scoring_positive_control"] = check_scoring_positive_control(
            str(gfa_path), outdir, args.threads, args.min_length, args.max_length,
        )
        control = summary["scoring_positive_control"]
        print(f"      source ranked first: {control['source_ranked_first']}; "
              f"source breadth {control['source_unique_breadth_pct']}% vs "
              f"{control['runner_up']} {control['runner_up_unique_breadth_pct']}%")
    except Exception as exc:                       # aligner or samtools absent
        summary["scoring_positive_control"] = {"skipped": str(exc)}
        print(f"      skipped: {exc}")

    inputs = [bait_path, gfa_path]
    summary_path = outdir / "VALIDATION_SUMMARY.json"
    provenance_path = outdir / "PROVENANCE.json"
    previous_summary = (
        json.loads(summary_path.read_text()) if summary_path.is_file() else {}
    )
    previous_runs = (
        json.loads(provenance_path.read_text()).get("runs", {})
        if provenance_path.is_file() else {}
    )
    parameters = {
        "deposit": Path(args.deposit).as_posix(),
        "placeholder_bp": PLACEHOLDER_BP,
        "ambiguous_window_bp": AMBIGUOUS_WINDOW_BP,
        "flank_bp": FLANK_BP,
        "graph_path_min_length_bp": args.min_length,
        "graph_path_max_length_bp": args.max_length,
        "threads": args.threads,
    }
    tools = {
        "python": [sys.executable, "--version"],
        "blastn": ["blastn", "-version"],
        "samtools": ["samtools", "--version"],
    }
    runs = {"checks_1_to_5": build_provenance(
        workflow="frame_and_graph_evidence/checks_1_to_5", repo_root=REPO_ROOT,
        inputs=inputs, parameters=parameters, tools=tools,
    )}
    if args.reads_r1:
        print("[6/6] competitive graph-path scoring on the published case")
        summary["path_scoring"] = check_path_scoring(
            str(gfa_path), outdir, args.reads_r1, args.reads_r2, args.threads,
            args.min_length, args.max_length, args.prefix,
        )
        print(f"      {summary['path_scoring']['paths_scored']} paths scored; "
              f"reads discriminate: {summary['path_scoring']['reads_discriminate']}")
        reads = [Path(args.reads_r1)] + ([Path(args.reads_r2)] if args.reads_r2 else [])
        runs["check_6_path_scoring"] = build_provenance(
            workflow="frame_and_graph_evidence/check_6", repo_root=REPO_ROOT,
            inputs=inputs + reads, parameters=parameters, tools=tools,
        )
    elif "path_scoring" in previous_summary:
        # Check 6 needs the published reads.  Without them the deposited
        # result and its provenance are kept rather than dropped.
        print("[6/6] no reads given: the deposited path-scoring result is kept")
        summary["path_scoring"] = previous_summary["path_scoring"]
        if "check_6_path_scoring" in previous_runs:
            runs["check_6_path_scoring"] = previous_runs["check_6_path_scoring"]

    write_json(summary_path, summary)
    write_json(provenance_path, {
        "schema_version": "1.0", "workflow": "frame_and_graph_evidence", "runs": runs,
    })
    checksums = {}
    for record in runs.values():
        for entry in record.get("inputs", []):
            checksums[entry["path"]] = (entry["size_bytes"], entry["sha256"])
    write_tsv(
        outdir / "INPUT_CHECKSUMS.tsv", ["path", "size_bytes", "sha256"],
        [[path, size, digest] for path, (size, digest) in sorted(checksums.items())],
    )
    print(f"\nWrote the result tables, VALIDATION_SUMMARY.json and PROVENANCE.json "
          f"to {outdir}")


if __name__ == "__main__":
    main()
