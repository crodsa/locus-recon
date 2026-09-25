#!/usr/bin/env python3
"""Deterministic benchmark for span-completeness and truncation measurement.

The completeness layer answers a geometric question: given a locus projected
onto a local assembly, how many query bases does the chosen contig fail to
supply, and is the missing sequence present on another contig?  That question
has an exact answer whenever the assembly is constructed rather than observed,
so it is benchmarked here against constructed truth instead of a biological
experiment.

Each case embeds a real locus sequence (H. pylori 26695 gyrB, the same bait
used in the low-copy validation) in synthetic contigs cut at known offsets, in
real genomic flanking context. The number of locus bases withheld from the
chosen contig is the ground truth. Measurement runs the production code path:
``blastn`` in ``BLAST_FMT_STEP5`` against the constructed assembly,
``rank_local_blast_hit``, ``collapse_overlapping_hits``,
``describe_full_query_span`` and ``find_span_continuation``, plus the
``SPAN_CLIP_THRESHOLDS`` reporting gate used by ``qc._classify_confidence``.
The tier columns classify each reported sequence with the identity and query
coverage the pipeline's validation search returns for it.

Cases cover: no loss with interior and flush contig boundaries (false-positive
controls), truncation at either query end immediately below and at the
reporting threshold, simultaneous truncation of both ends, truncation on a
reverse-orientation contig (where contig boundaries and query ends swap), and
three split-locus geometries -- overlapping continuation, butt-joined
continuation with no anchor, and a lost interval genuinely absent from the
assembly.

Usage:
    # Construct the cases from the source genome and measure them
    python validation/run_completeness_benchmark.py \
        --bait <fasta with gyrB record> \
        --genome <GenBank or FASTA of the source genome> \
        --outdir <output directory>

    # Re-measure the deposited constructed assemblies; no genome download
    python validation/run_completeness_benchmark.py \
        --from-deposit validation/completeness-benchmark \
        --outdir validation/completeness-benchmark
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from locus_recon.blast import (  # noqa: E402
    collapse_overlapping_hits,
    describe_full_query_span,
    find_span_continuation,
    rank_local_blast_hit,
)
from locus_recon.io import classify_result_disposition  # noqa: E402
from locus_recon.qc import (  # noqa: E402
    assess_allele_quality,
    profile_bait_database,
)
from locus_recon.utils import BLAST_FMT_STEP5, SPAN_CLIP_THRESHOLDS  # noqa: E402

FLANK = 500
DECOY = 800
PROFILE_CACHE: dict = {}


def _benchmark_remap(length: int) -> dict:
    """Uniform strong read support, so only span geometry can move the tier."""
    return {
        "mean_depth": 120.0,
        "min_depth": 80,
        "max_depth": 160,
        "breadth_pct": 100.0,
        "pct_bases_lt5": 0.0,
        "pct_bases_lt10": 0.0,
        "mapped_reads": 4000,
        "mapped_pct": 99.5,
        "per_base_available": True,
        "tier_uncertain_base_fraction": 0.0,
        "interior_uncertain_base_fraction": 0.0,
        "interior_metric_used": True,
        "positions_evaluated": length,
        "mean_base_quality": 36.0,
        "mean_mapping_quality": 59.0,
        "strand_balance_pct": 98.0,
    }


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_fasta(path: str) -> dict:
    records, name, chunks = {}, None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    records[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name:
        records[name] = "".join(chunks)
    return records


def load_genome(path: str) -> str:
    """Return the genome sequence from a GenBank or FASTA file."""
    if path.endswith((".fa", ".fasta", ".fna")):
        return next(iter(read_fasta(path).values())).upper()
    from Bio import SeqIO

    return str(next(SeqIO.parse(path, "genbank")).seq).upper()


def locate(genome: str, locus: str) -> tuple:
    """Find the locus in the genome; return (start, end, strand), 0-based."""
    idx = genome.find(locus)
    if idx >= 0:
        return idx, idx + len(locus), "+"
    idx = genome.find(revcomp(locus))
    if idx >= 0:
        return idx, idx + len(locus), "-"
    raise SystemExit("Locus sequence not found in the supplied genome.")


def build_cases(locus: str, up: str, down: str, decoy: str) -> list:
    """Construct every benchmark case with its constructed ground truth.

    Truth fields are in query (locus) coordinates:
      truth_lost_prefix_bp / truth_lost_suffix_bp -- locus bases the chosen
      contig cannot supply at the 5' / 3' query end;
      truth_clipped_start_bp / truth_clipped_end_bp -- the same loss expressed
      at the contig's own boundaries, which swap under reverse orientation.
    """
    n = len(locus)
    cases = []

    def add(name, contigs, prefix, suffix, clip_start, clip_end,
            continuation, rationale):
        cases.append({
            "case": name,
            "contigs": contigs,
            "truth_lost_prefix_bp": prefix,
            "truth_lost_suffix_bp": suffix,
            "truth_clipped_start_bp": clip_start,
            "truth_clipped_end_bp": clip_end,
            "truth_clipped_bp": clip_start + clip_end,
            "truth_continuation": continuation,
            "rationale": rationale,
        })

    add("complete_internal",
        {"ctg_complete": up + locus + down},
        0, 0, 0, 0, None,
        "Locus complete with flanking context on both sides: no loss.")

    add("flush_start_no_loss",
        {"ctg_flush_start": locus + down},
        0, 0, 0, 0, None,
        "Locus starts at contig base 1 with nothing missing: a contig "
        "boundary at the locus edge must not be read as truncation.")

    add("flush_end_no_loss",
        {"ctg_flush_end": up + locus},
        0, 0, 0, 0, None,
        "Locus ends at the last contig base with nothing missing.")

    add("truncated_start_9",
        {"ctg_trunc_start_9": locus[9:] + down},
        9, 0, 9, 0, None,
        "9 bp of the query 5' end withheld: one base below the reporting "
        "threshold, so measured but not reported as truncation.")

    add("truncated_start_10",
        {"ctg_trunc_start_10": locus[10:] + down},
        10, 0, 10, 0, None,
        "10 bp withheld: exactly at the reporting threshold.")

    add("truncated_end_9",
        {"ctg_trunc_end_9": up + locus[:n - 9]},
        0, 9, 0, 9, None,
        "9 bp of the query 3' end withheld: below the reporting threshold.")

    add("truncated_end_50",
        {"ctg_trunc_end_50": up + locus[:n - 50]},
        0, 50, 0, 50, None,
        "50 bp of the query 3' end withheld.")

    add("truncated_start_100",
        {"ctg_trunc_start_100": locus[100:] + down},
        100, 0, 100, 0, None,
        "100 bp of the query 5' end withheld.")

    add("truncated_end_100",
        {"ctg_trunc_end_100": up + locus[:n - 100]},
        0, 100, 0, 100, None,
        "100 bp of the query 3' end withheld.")

    add("truncated_both_30_40",
        {"ctg_trunc_both": locus[30:n - 40]},
        30, 40, 30, 40, None,
        "Both query ends withheld simultaneously (30 bp and 40 bp); the "
        "contig holds only an interior fragment of the locus.")

    add("reverse_truncated_suffix_60",
        {"ctg_reverse": revcomp(up + locus[:n - 60])},
        0, 60, 60, 0, None,
        "Reverse-orientation contig missing 60 bp of the query 3' end. The "
        "loss appears at the contig START, so orientation handling is what "
        "maps the clip back to the correct query end.")

    add("split_overlap_100",
        {"ctg_split_a": up + locus[:1800],
         "ctg_split_b": locus[1700:] + down},
        0, n - 1800, 0, n - 1800,
        {"sseqid": "ctg_split_b", "recovered_bp": n - 1800, "overlap_bp": 100},
        "Locus split across two contigs sharing 100 bp of query overlap: the "
        "lost interval is present in the assembly with sequence available to "
        "anchor a join.")

    add("split_no_overlap",
        {"ctg_split_c": up + locus[:1500],
         "ctg_split_d": locus[1500:] + down},
        0, n - 1500, 0, n - 1500,
        {"sseqid": "ctg_split_d", "recovered_bp": n - 1500, "overlap_bp": 0},
        "Butt-joined split with no shared query overlap: the lost interval is "
        "recoverable but no anchor exists, which is a different curation "
        "situation from an overlapping split.")

    add("split_absent_from_assembly",
        {"ctg_split_e": up + locus[:1500], "ctg_decoy": decoy},
        0, n - 1500, 0, n - 1500, None,
        "The lost interval is genuinely absent from the assembly: an "
        "unrelated contig must not be reported as a continuation.")

    return cases


VALIDATION_OUTFMT = "6 qseqid sseqid length pident qcovs bitscore qlen slen evalue"


def validate_like_pipeline(allele: str, bait_path: str, tools: dict) -> tuple:
    """Identity, query coverage and exact-match status as step 6 computes them.

    The pipeline searches the reconstructed allele against the bait set and
    reads ``qcovs``, the share of the *allele* the best bait hit covers, so an
    exact but truncated allele is covered 100%.  Its truncation is reported by
    the span layer, not by this number.
    """
    with tempfile.TemporaryDirectory(prefix="completeness-val-") as temporary:
        query = os.path.join(temporary, "allele.fasta")
        with open(query, "w") as handle:
            handle.write(f">allele\n{allele}\n")
        completed = subprocess.run(
            [tools["blastn"], "-query", query, "-subject", bait_path,
             "-outfmt", VALIDATION_OUTFMT],
            check=True, capture_output=True, text=True,
        )
    hits = [line.split("\t") for line in completed.stdout.splitlines()
            if len(line.split("\t")) >= 9]
    if not hits:
        return 0.0, 0.0, False
    hits.sort(
        key=lambda parts: (float(parts[5]), float(parts[4]),
                           float(parts[3]), int(parts[2])),
        reverse=True,
    )
    best = hits[0]
    identity, qcov = float(best[3]), float(best[4])
    exact = (identity == 100.0 and qcov == 100.0
             and int(best[2]) == int(best[6]) == int(best[7]))
    return identity, qcov, exact


def run_case(case: dict, locus: str, bait_path: str, outdir: str, tools: dict) -> dict:
    """Measure one case through the production span code path."""
    workdir = os.path.join(outdir, "cases", case["case"])
    os.makedirs(workdir, exist_ok=True)
    asm = os.path.join(workdir, "assembly.fasta")
    with open(asm, "w") as handle:
        for name, seq in case["contigs"].items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                handle.write(seq[i:i + 70] + "\n")

    hits_path = os.path.join(workdir, "local_hits.tsv")
    # The BLAST database is a build artefact of one BLAST version; only the
    # assembly and the hit table are kept with the case.
    with tempfile.TemporaryDirectory(prefix="completeness-db-") as temporary:
        db = os.path.join(temporary, "assembly_db")
        subprocess.run(
            [tools["makeblastdb"], "-in", asm, "-dbtype", "nucl", "-out", db],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        subprocess.run(
            [tools["blastn"], "-query", bait_path, "-db", db, "-out", hits_path,
             "-outfmt", BLAST_FMT_STEP5, "-evalue", "1e-10"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )

    hits = []
    with open(hits_path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 13:
                hits.append(rank_local_blast_hit(parts))
    if not hits:
        return {"case": case["case"], "measurement_ok": False,
                "note": "no BLAST hit in constructed assembly"}

    regions = collapse_overlapping_hits(hits)
    best = regions[0]
    span = describe_full_query_span(best, best["slen"])
    continuation = find_span_continuation(span, regions, best)

    # The reconstruction the pipeline would report is the projected span
    # extracted from the chosen contig.  Detection of a clip or of a
    # continuation must not alter it: this tool reports assembly geometry and
    # never joins contigs.
    contig_seq = case["contigs"][best["sseqid"]]
    extracted = contig_seq[span["start"]:span["end"]]
    expected_bp = best["qlen"] - case["truth_clipped_bp"]
    # Truth at base level: the part of the locus the chosen contig carries, in
    # query orientation.  The reported reconstruction must equal it exactly.
    oriented = revcomp(extracted) if best["sstart"] > best["send"] else extracted
    expected_sequence = locus[
        case["truth_lost_prefix_bp"]:len(locus) - case["truth_lost_suffix_bp"]
    ]
    identical_to_truth = oriented == expected_sequence

    # Effect of the completeness layer on the tier, in isolation: the same
    # sequence and read evidence classified with and without span metrics,
    # with identity and coverage computed as the pipeline's validation step
    # computes them.
    identity, qcov, exact_known = validate_like_pipeline(oriented, bait_path, tools)
    qc_with = assess_allele_quality(
        oriented, PROFILE_CACHE["profile"], identity, qcov, "benchmark",
        remap_metrics=_benchmark_remap(len(extracted)),
        ambiguity_metrics={"ambiguous": False, "bitscore_margin": 100.0},
        expect_cds=False,
        span_metrics=dict(span, continuation=continuation),
        exact_known_allele=exact_known,
    )
    qc_without = assess_allele_quality(
        oriented, PROFILE_CACHE["profile"], identity, qcov, "benchmark",
        remap_metrics=_benchmark_remap(len(extracted)),
        ambiguity_metrics={"ambiguous": False, "bitscore_margin": 100.0},
        expect_cds=False,
        exact_known_allele=exact_known,
    )
    span_flags = [
        f for f in qc_with["flags"]
        if f.startswith(("ALLELE_SPAN_CLIPPED_AT_CONTIG_END",
                         "LOCUS_SPLIT_ACROSS_CONTIGS"))
    ]

    unclaimed = span["unclaimed_query"]
    measured_prefix = next((hi - lo + 1 for lo, hi in unclaimed if lo == 1), 0)
    measured_suffix = next(
        (hi - lo + 1 for lo, hi in unclaimed if hi == best["qlen"]), 0
    )
    reported = span["clipped_bp"] >= SPAN_CLIP_THRESHOLDS["min_flag_bp"]

    row = {
        "case": case["case"],
        "best_contig": best["sseqid"],
        "orientation": "-" if best["sstart"] > best["send"] else "+",
        "aligned_bp": best["length"],
        "pident": round(best["pident"], 2),
        "truth_clipped_start_bp": case["truth_clipped_start_bp"],
        "measured_clipped_start_bp": span["clipped_start_bp"],
        "truth_clipped_end_bp": case["truth_clipped_end_bp"],
        "measured_clipped_end_bp": span["clipped_end_bp"],
        "truth_clipped_bp": case["truth_clipped_bp"],
        "measured_clipped_bp": span["clipped_bp"],
        "truth_lost_prefix_bp": case["truth_lost_prefix_bp"],
        "measured_lost_prefix_bp": measured_prefix,
        "truth_lost_suffix_bp": case["truth_lost_suffix_bp"],
        "measured_lost_suffix_bp": measured_suffix,
        "truncation_reported": reported,
        "truncation_reportable_by_truth":
            case["truth_clipped_bp"] >= SPAN_CLIP_THRESHOLDS["min_flag_bp"],
        "truth_continuation_contig":
            (case["truth_continuation"] or {}).get("sseqid", ""),
        "measured_continuation_contig":
            (continuation or {}).get("sseqid", ""),
        "truth_recovered_bp":
            (case["truth_continuation"] or {}).get("recovered_bp", 0),
        "measured_recovered_bp": (continuation or {}).get("recovered_bp", 0),
        "truth_overlap_bp":
            (case["truth_continuation"] or {}).get("overlap_bp", 0),
        "measured_overlap_bp": (continuation or {}).get("overlap_bp", 0),
        "span_flags": "; ".join(span_flags),
        "validation_identity_pct": identity,
        "validation_qcov_pct": qcov,
        "sequence_confidence_with_span": qc_with["sequence_confidence"],
        "sequence_confidence_without_span": qc_without["sequence_confidence"],
        "disposition_with_span":
            classify_result_disposition("SUCCESS", qc_with["confidence"]),
        "disposition_without_span":
            classify_result_disposition("SUCCESS", qc_without["confidence"]),
        "reconstructed_bp": len(extracted),
        "expected_reconstructed_bp": expected_bp,
        "reconstruction_identical_to_truth": identical_to_truth,
        "rationale": case["rationale"],
    }
    checks = [
        row["measured_clipped_start_bp"] == row["truth_clipped_start_bp"],
        row["measured_clipped_end_bp"] == row["truth_clipped_end_bp"],
        row["measured_lost_prefix_bp"] == row["truth_lost_prefix_bp"],
        row["measured_lost_suffix_bp"] == row["truth_lost_suffix_bp"],
        row["measured_continuation_contig"] == row["truth_continuation_contig"],
        row["measured_recovered_bp"] == row["truth_recovered_bp"],
        row["measured_overlap_bp"] == row["truth_overlap_bp"],
        row["truncation_reported"] == row["truncation_reportable_by_truth"],
        identical_to_truth,
        bool(span_flags) == row["truncation_reportable_by_truth"],
    ]
    row["exact"] = all(checks)
    row["measurement_ok"] = True
    return row


def cases_from_deposit(deposit: str) -> tuple:
    """Rebuild the case definitions from a deposited benchmark.

    The constructed assemblies are deposited with their construction truth, so
    the measurement can be repeated without the source genome.  Returns the
    locus sequence, its record name, the case list and the construction facts
    (strand and flank) that only the genome could otherwise supply.
    """
    deposited = json.load(
        open(os.path.join(deposit, "completeness_benchmark_results.json"))
    )
    query = read_fasta(os.path.join(deposit, "benchmark_query.fasta"))
    locus = next(iter(query.values())).upper()
    cases = []
    for record in deposited["cases"]:
        contigs = read_fasta(os.path.join(deposit, "cases", record["case"], "assembly.fasta"))
        continuation = None
        if record["truth_continuation_contig"]:
            continuation = {
                "sseqid": record["truth_continuation_contig"],
                "recovered_bp": record["truth_recovered_bp"],
                "overlap_bp": record["truth_overlap_bp"],
            }
        cases.append({
            "case": record["case"],
            "contigs": {name: seq.upper() for name, seq in contigs.items()},
            "truth_lost_prefix_bp": record["truth_lost_prefix_bp"],
            "truth_lost_suffix_bp": record["truth_lost_suffix_bp"],
            "truth_clipped_start_bp": record["truth_clipped_start_bp"],
            "truth_clipped_end_bp": record["truth_clipped_end_bp"],
            "truth_clipped_bp": record["truth_clipped_bp"],
            "truth_continuation": continuation,
            "rationale": record["rationale"],
        })
    summary = deposited["summary"]
    facts = {
        "locus_record": summary["locus_record"],
        "genome_locus_strand": summary["genome_locus_strand"],
        "flank_bp": summary["flank_bp"],
    }
    return locus, cases, facts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bait",
                        help="FASTA containing the locus record used as query")
    parser.add_argument("--locus-id", default="gyrB",
                        help="Substring identifying the record in --bait")
    parser.add_argument("--genome",
                        help="GenBank/FASTA genome supplying real flanks")
    parser.add_argument("--from-deposit", metavar="DIR",
                        help="re-measure the constructed assemblies deposited in DIR")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--blastn", default="blastn")
    args = parser.parse_args()
    if bool(args.from_deposit) == bool(args.bait or args.genome):
        parser.error("give either --from-deposit, or --bait and --genome")
    if not args.from_deposit and not (args.bait and args.genome):
        parser.error("--bait and --genome are required together")

    os.makedirs(args.outdir, exist_ok=True)
    if args.from_deposit:
        locus, cases, facts = cases_from_deposit(args.from_deposit)
        key = facts["locus_record"]
        record_name = key.split("|")[0]
    else:
        baits = read_fasta(args.bait)
        key = next((k for k in baits if args.locus_id in k), None)
        if key is None:
            raise SystemExit(f"No record matching {args.locus_id!r} in {args.bait}")
        locus = baits[key].upper()
        record_name = args.locus_id

    bait_path = os.path.join(args.outdir, "benchmark_query.fasta")
    with open(bait_path, "w") as handle:
        handle.write(f">{record_name}\n{locus}\n")

    PROFILE_CACHE["profile"] = profile_bait_database(bait_path)

    if args.from_deposit:
        strand, flank = facts["genome_locus_strand"], facts["flank_bp"]
    else:
        genome = load_genome(args.genome)
        start, end, strand = locate(genome, locus)
        if strand == "+":
            up, down = genome[start - FLANK:start], genome[end:end + FLANK]
        else:
            up = revcomp(genome[end:end + FLANK])
            down = revcomp(genome[start - FLANK:start])
        decoy_start = (end + 200_000) % (len(genome) - DECOY)
        decoy = genome[decoy_start:decoy_start + DECOY]
        cases = build_cases(locus, up, down, decoy)
        flank = FLANK

    tools = {"makeblastdb": args.makeblastdb, "blastn": args.blastn}
    rows = [run_case(case, locus, bait_path, args.outdir, tools) for case in cases]

    columns = list(rows[0].keys())
    tsv = os.path.join(args.outdir, "completeness_benchmark_results.tsv")
    with open(tsv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    exact = sum(1 for r in rows if r.get("exact"))
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"),
        "locus_record": key,
        "locus_length_bp": len(locus),
        "genome_locus_strand": strand,
        "flank_bp": flank,
        "reporting_threshold_bp": SPAN_CLIP_THRESHOLDS["min_flag_bp"],
        "cases": len(rows),
        "cases_exact": exact,
        "cases_measured": sum(1 for r in rows if r.get("measurement_ok")),
        "max_abs_error_bp": max(
            abs(r.get("measured_clipped_bp", 0) - r.get("truth_clipped_bp", 0))
            for r in rows
        ),
        "false_positive_cases": [
            r["case"] for r in rows
            if r.get("truth_clipped_bp") == 0 and r.get("measured_clipped_bp")
        ],
        "false_negative_cases": [
            r["case"] for r in rows
            if r.get("truth_clipped_bp") and not r.get("measured_clipped_bp")
        ],
        "cases_reconstruction_identical_to_truth": sum(
            1 for r in rows if r.get("reconstruction_identical_to_truth")),
        "tier_changed_by_span_metrics": [
            r["case"] for r in rows
            if r.get("sequence_confidence_with_span")
            != r.get("sequence_confidence_without_span")
        ],
        "blastn_outfmt": BLAST_FMT_STEP5,
    }
    with open(os.path.join(args.outdir, "completeness_benchmark_results.json"),
              "w") as handle:
        json.dump({"summary": summary, "cases": rows}, handle, indent=2)
        handle.write("\n")

    print(f"{exact}/{len(rows)} cases exact; TSV: {tsv}")
    return 0 if exact == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
