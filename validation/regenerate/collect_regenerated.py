#!/usr/bin/env python3
"""Assemble the package of a regeneration of the validation analyses.

``regenerate_validation.sh`` reruns the analyses and leaves their outputs in a
working directory.  This script turns those outputs into the package that is
kept, in four commands:

``evidence``  copies, for every Locus-Recon run, the files that let a result be
              checked: the candidate FASTA, the QC report, the per-base table
              ``allele_remap.support.tsv``, the reports, logs, parameters and
              the local graph.  Alignments, reads, indexes and any file above
              5 MB stay in the working directory and are listed.
``records``   writes ``records/`` in the layout of ``validation/`` in the
              repository, so that copying it over ``validation/`` replaces the
              records of the earlier release.
``compare``   writes ``old_vs_new.tsv``: every value of every record of the
              earlier release beside the regenerated value.
``manifest``  writes ``PROVENANCE.json`` and ``SHA256SUMS.txt``.

Only the standard library is used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path

MAX_BYTES = 5 * 1024 * 1024
STRAINS = ("Hpfe0001", "Hpfe0002", "Hpfe0003", "Hpfe0004", "Hpfe0006")
RUN_PATTERNS = (
    "locus_recon_report_*.tsv", "run_manifest.json", "locus_recon_batch_*.log",
    "*_accepted_alleles.fasta", "*_hold_candidates.fasta",
    "*_review_required_candidates.fasta", "*_reconstructed_alleles.fasta",
    "*_reconstructed_candidates.fasta",
)
SAMPLE_PATTERNS = (
    "*_reconstructed.fasta", "*_qc_report.txt", "*_recon.log",
    "allele_remap.support.tsv", "allele_remap.flagstat.txt", "allele_remap_bwa.log",
    "allele_remap_singleton_bwa.log", "allele_region.bed", "target_regions.bed",
    "blast_step1.tsv", "blast_step5.tsv", "blast_validation.tsv", "bwa_pipe.log",
    "mapping_flagstat.txt",
)
SPADES_FILES = (
    "assembly_graph_after_simplification.gfa", "contigs.fasta", "scaffolds.fasta",
    "params.txt", "spades.log",
)
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ------------------------------------------------------------------ helpers

def digest_file(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path, delimiter=None):
    path = Path(path)
    delimiter = delimiter or ("," if path.suffix == ".csv" else "\t")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_table(path, header, rows, lineterminator="\n"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator=lineterminator)
        writer.writerow(header)
        writer.writerows(rows)


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def read_fasta(path):
    records, name, parts = [], None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                name, parts = line[1:], []
            elif line:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def copy_file(source, destination, skipped):
    source, destination = Path(source), Path(destination)
    if not source.is_file():
        return False
    size = source.stat().st_size
    if size > MAX_BYTES:
        skipped.append((str(source), size, "larger than 5 MB"))
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    return True


def copy_run(run_dir, destination, skipped):
    """Evidence files of one output directory of ``locus-recon``."""
    run_dir, destination = Path(run_dir), Path(destination)
    if not run_dir.is_dir():
        return 0
    copied = 0
    for pattern in RUN_PATTERNS:
        for path in sorted(run_dir.glob(pattern)):
            copied += copy_file(path, destination / path.name, skipped)
    for sample_dir in sorted(path for path in run_dir.iterdir() if path.is_dir()):
        if sample_dir.name in ("bait_database", "validation"):
            continue
        for pattern in SAMPLE_PATTERNS:
            for path in sorted(sample_dir.glob(pattern)):
                copied += copy_file(path, destination / sample_dir.name / path.name, skipped)
        for name in SPADES_FILES:
            copied += copy_file(sample_dir / "spades_local" / name,
                                destination / sample_dir.name / "spades_local" / name, skipped)
    return copied


def assembly_statistics(path):
    lengths = sorted((len(sequence) for _, sequence in read_fasta(path)), reverse=True)
    total, running, n50 = sum(lengths), 0, 0
    for length in lengths:
        running += length
        if running * 2 >= total:
            n50 = length
            break
    return len(lengths), total, n50


def canonical_sequences(path):
    return sorted(min(sequence.upper(), sequence.upper().translate(COMPLEMENT)[::-1])
                  for _, sequence in read_fasta(path))


# ----------------------------------------------------------------- evidence

def command_evidence(args):
    work, out = args.work, args.package / "evidence"
    if out.exists():
        shutil.rmtree(out)
    skipped = []

    # simulated benchmark: the run that is the record, and its repetition
    mock = work / "mock" / "mock_results"
    copy_run(mock, out / "mock", skipped)
    for name in ("validation_results.tsv", "validation_summary.json",
                 "workflow_provenance.json"):
        copy_file(mock / "validation" / name, out / "mock" / "validation" / name, skipped)
    repeat = work / "mock" / "mock_results_repeat"
    rows = []
    if repeat.is_dir():
        for path in sorted((out / "mock").rglob("*")):
            if not path.is_file():
                continue
            relative = path.relative_to(out / "mock")
            other = repeat / relative
            first = digest_file(path)
            second = digest_file(other) if other.is_file() else ""
            rows.append([relative.as_posix(), path.stat().st_size, first, second,
                         "yes" if first == second else "no"])
        write_table(out / "mock_run_and_repeat_checksums.tsv",
                    ["file", "size_bytes", "sha256_run", "sha256_repeat", "identical"], rows)
        for name in ("locus_recon_report_mockLocus.tsv", "run_manifest.json"):
            copy_file(repeat / name, out / "mock_repeat" / name, skipped)
        for name in ("validation_results.tsv", "validation_summary.json"):
            copy_file(repeat / "validation" / name, out / "mock_repeat" / "validation" / name,
                      skipped)

    # completeness benchmark measured again on the deposited assemblies
    for name in ("completeness_benchmark_results.tsv", "completeness_benchmark_results.json"):
        copy_file(work / "completeness" / "from_deposit" / name,
                  out / "completeness_from_deposited_assemblies" / name, skipped)

    # H. pylori: every case of stages B, C and D of the tier calibration
    drafts = []
    for strain in STRAINS:
        strain_dir = work / "tier" / strain
        if not strain_dir.is_dir():
            continue
        for report in sorted(strain_dir.glob("*/locus_recon_report_*.tsv")):
            case = report.parent.name
            copy_run(report.parent, out / "hpylori" / strain / case, skipped)
            copy_file(strain_dir / f"samples_{case}.tsv",
                      out / "hpylori" / strain / case / f"samples_{case}.tsv", skipped)
        for draft in sorted(strain_dir.glob("spades*")):
            contigs = draft / "contigs.fasta"
            if not contigs.is_file():
                continue
            copy_file(draft / "params.txt",
                      out / "hpylori" / strain / "drafts" / draft.name / "params.txt", skipped)
            count, total, n50 = assembly_statistics(contigs)
            row = [strain, draft.name, count, total, n50, digest_file(contigs),
                   digest_file(draft / "assembly_graph_after_simplification.gfa")
                   if (draft / "assembly_graph_after_simplification.gfa").is_file() else ""]
            earlier = (args.earlier_tier_work / strain / draft.name / "contigs.fasta"
                       if args.earlier_tier_work else None)
            if earlier is not None and earlier.is_file():
                earlier_count, earlier_total, earlier_n50 = assembly_statistics(earlier)
                row += [earlier_count, earlier_total, earlier_n50, digest_file(earlier),
                        "yes" if digest_file(earlier) == row[5] else "no",
                        "yes" if canonical_sequences(earlier) == canonical_sequences(contigs)
                        else "no"]
            else:
                row += ["", "", "", "", "", ""]
            drafts.append(row)
    if drafts:
        write_table(out / "hpylori" / "draft_assemblies.tsv", [
            "strain", "draft", "contigs", "total_bp", "n50_bp", "contigs_fasta_sha256",
            "graph_gfa_sha256", "earlier_contigs", "earlier_total_bp", "earlier_n50_bp",
            "earlier_contigs_fasta_sha256", "byte_identical_to_earlier",
            "same_contig_sequences_as_earlier"], drafts)

    # LIBA-6656: the two runs of the standard workflow
    liba = work / "liba6656_gfa_rerun"
    for run in sorted(liba.glob("run_*_threads")):
        copy_run(run / "output", out / "liba6656" / run.name, skipped)
        copy_file(run / "samples.tsv", out / "liba6656" / run.name / "samples.tsv", skipped)
    for name in ("tcdB_terminal_paths_8000_9000.tsv", "RELEASE_VERIFICATION.json"):
        copy_file(liba / name, out / "liba6656" / name, skipped)
    paths = liba / "tcdB_terminal_paths_8000_9000.fasta"
    if paths.is_file():
        (out / "liba6656").mkdir(parents=True, exist_ok=True)
        (out / "liba6656" / "tcdB_terminal_paths_8000_9000.fasta.sha256").write_text(
            f"{digest_file(paths)}  tcdB_terminal_paths_8000_9000.fasta\n")

    # copy-number workflows: results, provenance and logs (no alignments)
    def copy_small(directory, target, patterns):
        for pattern in patterns:
            for path in sorted(Path(directory).glob(pattern)):
                copy_file(path, out / "copy_number" / target / path.name, skipped)

    copy_small(work / "aphA1" / "results", "aphA1", ("*.tsv", "*.json", "*.log"))
    copy_small(work / "low_copy_23S" / "results", "low_copy_23S", ("*.tsv", "*.json"))
    copy_small(work / "low_copy_23S" / "fixed_baits", "low_copy_23S/fixed_baits", ("*.fasta",))
    for strain in STRAINS:
        copy_small(work / "low_copy_23S" / "samples" / strain, f"low_copy_23S/samples/{strain}",
                   ("fixed_baits_vs_draft.tsv", "*.mapping.log", "spades.command.txt"))
    copy_small(work / "constructed_copy_series" / "results", "constructed_copy_series",
               ("*.tsv", "*.json"))
    copy_small(work / "liba6656_depth" / "results", "liba6656_depth",
               ("*.tsv", "*.json", "*.log"))
    copy_small(work / "overlap_null", "overlap_null", ("*.csv", "*.json", "*.tsv", "*.log"))

    copy_file(work / "frame_graph" / "thread_count_comparison.tsv",
              out / "frame_graph" / "thread_count_comparison.tsv", skipped)

    exactness_check(args, out)
    write_table(args.package / "FILES_LEFT_ON_HOST.tsv", ["path", "size_bytes", "reason"],
                [[path, size, reason] for path, size, reason in sorted(skipped)])
    files = [path for path in out.rglob("*") if path.is_file()]
    print(f"evidence: {len(files)} files, {sum(path.stat().st_size for path in files)} bytes; "
          f"{len(skipped)} file(s) above 5 MB left on the host")


def exactness_check(args, evidence):
    """Compare every kept candidate with its truth sequence, base for base.

    The check does not use the runners: a candidate is exact when it is the
    truth sequence or its reverse complement.  The label that the runner
    recorded for the same case is written beside the verdict.
    """
    import gzip

    def reverse(sequence):
        return sequence.translate(COMPLEMENT)[::-1]

    def verdict(candidate, truths):
        candidate = candidate.upper()
        for name, truth in truths:
            truth = truth.upper()
            if candidate in (truth, reverse(truth)):
                return "yes", name, len(truth)
        if not truths:
            return "no truth sequence", "", ""
        return "no", truths[0][0], len(truths[0][1])

    rows = []
    repo, work = args.repo, args.work
    if repo is None:
        return
    # simulated benchmark
    mock = repo / "validation" / "mock_dataset"
    alleles = {header.split()[0]: sequence
               for header, sequence in read_fasta(mock / "truth" / "truth_alleles.fasta")}
    manifest = json.loads((mock / "dataset_manifest.json").read_text())
    recorded = {}
    table = evidence / "mock" / "validation" / "validation_results.tsv"
    if table.is_file():
        recorded = {row["sample_id"]: row["exact_truth_match"] for row in read_table(table)}
    for entry in manifest["samples"]:
        sample = entry["sample_id"]
        candidates = sorted((evidence / "mock" / sample).glob("*_reconstructed.fasta"))
        truths = [(entry["truth_allele"], alleles[entry["truth_allele"]])] \
            if entry.get("truth_allele") else []
        label = {"True": "yes", "False": "no"}.get(recorded.get(sample, ""), "")
        if not candidates:
            rows.append(["simulated benchmark", sample, "", 0, "", "truth/truth_alleles.fasta",
                         entry.get("truth_allele") or "", len(truths[0][1]) if truths else "",
                         "no candidate", label, "yes" if label in ("no", "") else "no"])
            continue
        sequence = read_fasta(candidates[0])[0][1]
        exact, name, length = verdict(sequence, truths)
        rows.append(["simulated benchmark", sample,
                     candidates[0].relative_to(evidence.parent).as_posix(), len(sequence),
                     hashlib.sha256(sequence.upper().encode()).hexdigest(),
                     "truth/truth_alleles.fasta", name, length, exact, label,
                     "yes" if exact == label else "no"])

    # H. pylori: annotated genes of the closed genomes
    def annotation(accession, kind):
        found = sorted((work / "tier" / "downloads").glob(
            f"{accession}_*_{kind}_from_genomic.fna.gz"))
        records, header, parts = [], None, []
        if not found:
            return records
        with gzip.open(found[0], "rt") as handle:
            for line in handle:
                if line.startswith(">"):
                    if header:
                        records.append((header, "".join(parts)))
                    header, parts = line[1:].strip(), []
                else:
                    parts.append(line.strip())
        if header:
            records.append((header, "".join(parts)))
        return records

    labels = {}
    tier = args.driver_records / "tier-calibration"
    for name in ("closed_genomes", "divergent_loci", "depth_series", "depth_series_isolate"):
        path = tier / f"tier_calibration_{name}.tsv"
        if path.is_file():
            for row in read_table(path):
                labels[row["case"]] = row["matches_truth"]
    manifest_path = repo / "validation" / "tier-calibration" / "closed_genome_manifest.tsv"
    for entry in read_table(manifest_path) if manifest_path.is_file() else []:
        strain, accession = entry["strain"], entry["assembly_accession"]
        cds, rna = annotation(accession, "cds"), annotation(accession, "rna")
        for case_dir in sorted((evidence / "hpylori" / strain).glob("*")):
            candidates = sorted(case_dir.glob(f"{strain}/*_reconstructed.fasta"))
            if case_dir.name == "drafts" or not case_dir.is_dir():
                continue
            locus, _, arm = case_dir.name.partition("_depth")
            if locus == "23S":
                truths = [(header.split()[0], sequence) for header, sequence in rna
                          if "[product=23S ribosomal RNA]" in header]
            else:
                truths = [(header.split()[0], sequence) for header, sequence in cds
                          if f"[gene={locus}]" in header]
            case = f"{strain} {locus}" + (f" {arm}" if arm else "")
            label = labels.get(case, "")
            source = f"{accession} {'rna' if locus == '23S' else 'cds'}_from_genomic.fna.gz"
            if not candidates:
                rows.append(["H. pylori closed genomes", case, "", 0, "", source, "", "",
                             "no candidate", label, "yes" if label in ("no", "") else "no"])
                continue
            sequence = read_fasta(candidates[0])[0][1]
            exact, name, length = verdict(sequence, truths)
            rows.append(["H. pylori closed genomes", case,
                         candidates[0].relative_to(evidence.parent).as_posix(), len(sequence),
                         hashlib.sha256(sequence.upper().encode()).hexdigest(), source, name,
                         length, exact, label, "yes" if exact == label else "no"])
    write_table(args.package / "exactness_check.tsv", [
        "analysis", "case", "candidate_file", "candidate_bp", "candidate_sha256",
        "truth_source", "truth_record", "truth_bp", "identical_to_truth",
        "label_recorded_by_the_runner", "check_agrees_with_label"], rows)
    print(f"exactness check: {len(rows)} cases, "
          f"{sum(row[8] == 'yes' for row in rows)} identical to the truth, "
          f"{sum(row[10] == 'no' for row in rows)} disagreement(s) with the recorded label")


# ------------------------------------------------------------------ records

def relative_report(report, run_dir, destination):
    """Copy a batch report with paths written relative to its directory."""
    run_dir = Path(run_dir).resolve()
    prefix = str(run_dir) + "/"
    with open(report, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    for row in rows[1:]:
        for index, value in enumerate(row):
            if value.startswith(prefix):
                row[index] = value[len(prefix):]
    Path(destination).parent.mkdir(parents=True, exist_ok=True)
    with open(destination, "w", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


def compact(provenance, **extra):
    software = provenance["software"]
    return {
        "schema_version": "1.0", "workflow": provenance["workflow"],
        "created_utc": provenance["created_utc"],
        "software": {key: software[key] for key in ("name", "version", "source_sha256")},
        "parameters": provenance["parameters"],
        "tools": {name: record["version"] for name, record in provenance["tools"].items()},
        **extra,
    }


def mock_reproduction(records, work, package):
    """BENCHMARK_REPRODUCTION.md of the simulated benchmark, from the run itself."""
    target = records / "reference_results"
    report = read_table(target / "locus_recon_report_mockLocus.tsv")
    checked = {row["sample_id"]: row for row in read_table(target / "validation_results.tsv")}
    summary = json.loads((target / "validation_summary.json").read_text())
    provenance = json.loads(
        (work / "mock" / "mock_results" / "validation" / "workflow_provenance.json").read_text())
    by_status = {}
    for row in report:
        by_status[row["workflow_status"]] = by_status.get(row["workflow_status"], 0) + 1
    tools = {name: record["version"] for name, record in provenance["tools"].items()}
    timed = [float(row["elapsed_sec"]) for row in report
             if row["workflow_status"] == "SUCCESS" and row["elapsed_sec"]]
    lines = [
        "# Benchmark reproduction record", "",
        "The compact files in this directory record the behaviour of the bundled synthetic",
        "benchmark. They are reference results for this deterministic panel, not universal",
        "estimates of biological performance.", "",
        f"The run recorded here was made from a clean results directory with Locus-Recon "
        f"{provenance['software']['version']} on {provenance['created_utc'][:10]}.", "",
        "```bash", "python validation/run_mock_validation.py \\",
        f"  --force \\", f"  --threads {provenance['parameters']['threads']} \\",
        f"  --memory-gb {provenance['parameters']['memory_gb']}", "```", "",
        "Environment:", "",
    ]
    lines += [f"- {name}: {version}" for name, version in sorted(tools.items())]
    lines += [
        "- package source digest (`source_sha256`, see `validation/common.py`):",
        f"  `{provenance['software']['source_sha256']}`", "", "Outcome:", "",
        f"- {len(report)} samples evaluated;",
        "- workflow status: " + ", ".join(
            f"{count} `{status}`" for status, count in sorted(by_status.items())) + ";",
        f"- core criteria passed: {sum(summary['criteria'].values())} of "
        f"{len(summary['criteria'])} (`all_core_criteria_passed`: "
        f"{str(summary['all_core_criteria_passed']).lower()});",
    ]
    for row in report:
        check = checked[row["sample_id"]]
        lines.append(
            f"- `{row['sample_id']}`: status `{row['workflow_status']}`, tier "
            f"`{row['sequence_confidence'] or 'none'}`, disposition "
            f"`{row['result_disposition'] or 'none'}`, identical to the truth "
            f"{check['exact_truth_match'].lower()}, remap depth {check['remap_mean_depth']}x, "
            f"mixed sites {check['mixed_site_count']} (median alternative fraction "
            f"{check['median_mixed_fraction']}), strand-biased sites "
            f"{check['strand_biased_sites']}, ambiguity flagged "
            f"{check['ambiguity_flagged'].lower()};")
    if timed:
        lines += ["", f"Timing: the {len(timed)} samples that entered read recruitment took a "
                  f"median of {statistics.median(timed):.2f} s each (range {min(timed):.1f} to "
                  f"{max(timed):.1f} s), {sum(timed):.1f} s in total. These figures describe "
                  "this run on this machine."]
    table = package / "evidence" / "mock_run_and_repeat_checksums.tsv"
    if table.is_file():
        rows = read_table(table)
        same = sum(row["identical"] == "yes" for row in rows)
        lines += ["", f"A second run of the same command was compared file by file with this "
                  f"one: {same} of {len(rows)} kept files were byte-identical "
                  "(`evidence/mock_run_and_repeat_checksums.tsv` in the regeneration package "
                  "lists every file)."]
    lines += ["", "In `locus_recon_report_mockLocus.tsv`, `allele_file` is shown relative to "
              "the results directory; the pipeline itself writes absolute paths.", ""]
    (target / "BENCHMARK_REPRODUCTION.md").write_text("\n".join(lines))


def command_records(args):
    work, records, repo, old = args.work, args.package / "records", args.repo, args.old
    records.mkdir(parents=True, exist_ok=True)
    built = []

    # simulated benchmark and completeness: written by the driver; add the text record
    for name in ("reference_results", "completeness-benchmark", "constructed-copy-series"):
        source = args.driver_records / name
        if source.is_dir():
            if (records / name).exists():
                shutil.rmtree(records / name)
            shutil.copytree(source, records / name)
            built.append(name)
    if (records / "reference_results").is_dir():
        mock_reproduction(records, work, args.package)

    # tier calibration
    source = args.driver_records / "tier-calibration"
    if source.is_dir():
        if (records / "tier-calibration").exists():
            shutil.rmtree(records / "tier-calibration")
        shutil.copytree(source, records / "tier-calibration")
        built.append("tier-calibration")

    # aphA1
    results = work / "aphA1" / "results"
    if (results / "workflow_provenance.json").is_file():
        target = records / "aphA1-copy-number"
        target.mkdir(parents=True, exist_ok=True)
        # The record of 1.0.0 shows the four numeric columns with three
        # decimals and Unix line ends, and evaluation.json in the key order
        # of that record; the files of the runner are kept under evidence/.
        rows = read_table(results / "aphA1_copy_number.tsv")
        for row in rows:
            for column in ("dosage_estimate", "ratio_ci_low", "ratio_ci_high",
                           "ambiguity_index"):
                if row.get(column) not in (None, ""):
                    row[column] = f"{float(row[column]):.3f}"
        write_table(target / "results.tsv", list(rows[0]), [list(row.values()) for row in rows])
        provenance = json.loads((results / "workflow_provenance.json").read_text())
        evaluation = json.loads((results / "evaluation.json").read_text())
        earlier_evaluation = json.loads(
            (old / "aphA1-copy-number" / "evaluation.json").read_text())
        order = [key for key in earlier_evaluation if key in evaluation]
        order += [key for key in evaluation if key not in order]
        (target / "evaluation.json").write_text(
            json.dumps({key: evaluation[key] for key in order}, indent=2) + "\n")
        for name in ("results.tsv", "evaluation.json"):
            same = (target / name).read_bytes() == (old / "aphA1-copy-number" / name).read_bytes()
            print(f"aphA1-copy-number/{name}: same bytes as the 1.0.0 record: {same}")
        full = json.loads((results / "aphA1_depth_ratio.json").read_text())
        parameters = provenance["parameters"]
        write_json(target / "run_summary.json", {
            "schema_version": "1.0", "workflow": provenance["workflow"],
            "created_utc": provenance["created_utc"],
            "all_prespecified_criteria_passed": evaluation["all_prespecified_criteria_passed"],
            "locus_region": parameters["aphA1_region"], "reference": "GCA_019458485.1",
            "reference_sha256": next(entry["sha256"] for entry in provenance["inputs"]
                                     if entry["path"].endswith(".fasta")),
            "runs": sorted(row["run"] for row in rows),
            "parameters": {
                "bootstrap_block_size": full[0]["result"].get("bootstrap_block_size"),
                "bootstrap_resamples": parameters["n_boot"],
                "call_threshold": parameters["call_threshold"],
                "max_backbone_bootstrap_blocks": parameters["max_backbone_bootstrap_blocks"],
                "min_bq": parameters["min_bq"], "min_mq": parameters["min_mq"],
                "seed": parameters["seed"], "threads": parameters["threads"],
            },
            "software": {key: provenance["software"][key]
                         for key in ("name", "version", "source_sha256")},
            "tools": {
                "bwa-mem2": provenance["tools"]["aligner"]["version"],
                "python": provenance["tools"]["python"]["version"].replace("Python ", ""),
                "samtools": provenance["tools"]["samtools"]["version"].replace("samtools ", ""),
            },
        })
        built.append("aphA1-copy-number")

    # low-copy 23S
    results = work / "low_copy_23S" / "results"
    if (results / "low_copy_results.tsv").is_file():
        target = records / "low-copy-23S"
        target.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(results / "low_copy_results.tsv", target / "results.tsv")
        for name in ("closed_truth_audit.json", "evaluation.json", "workflow_provenance.json"):
            if (results / name).is_file():
                shutil.copyfile(results / name, target / name)
        shutil.copyfile(repo / "validation" / "low-copy-23S" / "PRESPECIFICATION.md",
                        target / "PRESPECIFICATION.md")
        rows = read_table(target / "results.tsv")
        rrna = [row for row in rows if row["locus"] == "23S"]
        controls = [row for row in rows if row["locus"] == "gyrB"]
        earlier = json.loads((old / "low-copy-23S" / "run_summary.json").read_text())
        legacy = legacy_fields_unchanged(old / "low-copy-23S" / "results.tsv", rows, repo)
        if (target / "evaluation.json").is_file() and (target / "workflow_provenance.json").is_file():
            evaluation = json.loads((target / "evaluation.json").read_text())
            provenance = json.loads((target / "workflow_provenance.json").read_text())
            multicopy = sum(row["depth_call"] == "MULTICOPY_DEPTH" for row in rrna)
            false_positive = sum(row["depth_call"] == "MULTICOPY_DEPTH" for row in controls)
            panel = {
                "samples": evaluation["samples"],
                "two_copy_23S_multicopy_depth_calls": multicopy,
                "two_copy_sensitivity": evaluation["two_copy_sensitivity"],
                "single_copy_gyrB_false_positive_calls": false_positive,
                "legacy_depth_fields_unchanged": legacy["unchanged"],
                "prespecified_depth_criteria_passed": (
                    evaluation["criterion_two_copy_sensitivity_ge_0_80"]
                    and evaluation["criterion_zero_control_false_positives"]),
                "exact_two_copy_23S_consensus_calls": evaluation["exact_two_copy_consensus_calls"],
                "two_copy_23S_graph_lower_bounds": sum(
                    row["copy_number_kind"] == "LOWER_BOUND" for row in rrna),
                "exact_single_copy_gyrB_consensus_calls":
                    evaluation["exact_single_copy_control_consensus_calls"],
                "all_primary_criteria_passed": evaluation["all_primary_criteria_passed"],
            }
            ratios = [float(row["ratio"]) for row in rrna]
            review = [row["strain"] for row in rrna if row["copy_number_status"] == "REVIEW"]
            interpretation = (
                f"The aggregate-depth rule classified {multicopy} of {len(rrna)} two-copy 23S "
                f"loci as MULTICOPY_DEPTH and {false_positive} of {len(controls)} single-copy "
                f"gyrB controls as MULTICOPY_DEPTH. Observed 23S depth ratios were "
                f"{min(ratios):.3f}-{max(ratios):.3f}. The copy-number consensus was an exact, "
                f"supported count of two for {panel['exact_two_copy_23S_consensus_calls']} of "
                f"{len(rrna)} 23S loci, a lower bound for "
                f"{panel['two_copy_23S_graph_lower_bounds']}, and an exact, supported count of "
                f"one for {panel['exact_single_copy_gyrB_consensus_calls']} of {len(controls)} "
                f"gyrB controls. 23S loci with copy_number_status REVIEW: "
                f"{', '.join(review) or 'none'}. all_primary_criteria_passed is "
                f"{str(panel['all_primary_criteria_passed']).lower()}. The panel informed the "
                "design of the graph method; it is a biological demonstration, not an "
                "independent performance estimate."
            )
            write_json(target / "run_summary.json", {
                "schema_version": "1.0", "workflow": provenance["workflow"],
                "completed_utc": provenance["created_utc"],
                "software": {key: provenance["software"][key]
                             for key in ("name", "version", "source_sha256")},
                "panel": panel, "interpretation": interpretation,
                "legacy_depth_field_differences": legacy["differences"],
                "execution_note": (
                    "The five SPAdes --isolate drafts were assembled by the closed-genome "
                    "stage of the tier calibration of this regeneration (4 threads, 14 GB "
                    "memory cap) from the same checksum-verified read files and reused here, "
                    "as in the record of release 1.0.0: the workflow keeps an existing draft, "
                    "so memory_gb below applied to no assembly in this run. The GFA checksums "
                    "are recorded in workflow_provenance.json."),
                "panel_of_release_1_0_0": earlier.get("panel"),
                "parameters": {key: value for key, value in provenance["parameters"].items()
                               if key != "selection"},
                "result_sha256": digest_file(target / "results.tsv"),
                "evaluation_sha256": digest_file(target / "evaluation.json"),
                "truth_audit_sha256": digest_file(target / "closed_truth_audit.json"),
                "provenance_sha256": digest_file(target / "workflow_provenance.json"),
            })
        else:
            write_json(target / "run_summary.json", {
                "schema_version": "1.0", "workflow": "external-two-copy-23S",
                "complete": False,
                "note": "The workflow stopped at its check of the archived depth fields; "
                        "results.tsv holds the values it had computed.",
                "legacy_depth_fields_unchanged": legacy["unchanged"],
                "legacy_depth_field_differences": legacy["differences"],
            })
        built.append("low-copy-23S")

    # LIBA-6656 depth
    results = work / "liba6656_depth" / "results"
    if (results / "workflow_provenance.json").is_file():
        target = records / "liba6656-depth"
        target.mkdir(parents=True, exist_ok=True)
        for name in ("LIBA6656_tcdB_depth_result.tsv", "evaluation.json",
                     "workflow_provenance.json"):
            shutil.copyfile(results / name, target / name)
        row = read_table(results / "LIBA6656_tcdB_depth_result.tsv")[0]
        earlier_header = list(read_table(old / "liba6656-depth" / "result_summary.tsv")[0])
        write_table(target / "result_summary.tsv", earlier_header,
                    [[row.get(key, "") for key in earlier_header]], lineterminator="\r\n")
        provenance = json.loads((results / "workflow_provenance.json").read_text())
        names = {"LIBA6656_ST154.assembly.fasta", "ERR467623_1.fastq.gz",
                 "ERR467623_2.fastq.gz", "liba6656_depth_regions.tsv"}
        write_json(target / "run_summary.json", compact(
            provenance, complete=provenance.get("complete"),
            input_sha256={Path(entry["path"]).name: entry["sha256"]
                          for entry in provenance["inputs"]
                          if Path(entry["path"]).name in names},
            note="workflow_provenance.json is the record validation/workflows/"
                 "run_liba6656_depth.py wrote for this run: every input with its checksum, "
                 "the parameters and the tool versions. LIBA6656_tcdB_depth_result.tsv is "
                 "the full result table it lists; result_summary.tsv and evaluation.json "
                 "carry the reported fields."))
        built.append("liba6656-depth")

    # overlap null
    source = work / "overlap_null"
    if (source / "null_dose_overlap_control.json").is_file():
        target = records / "depth-copy-number" / "overlap_null"
        target.mkdir(parents=True, exist_ok=True)
        for name in ("null_dose_overlap_control.csv", "null_dose_overlap_control.json"):
            shutil.copyfile(source / name, target / name)
        built.append("depth-copy-number/overlap_null")

    # LIBA-6656 rerun of the standard workflow
    liba = work / "liba6656_gfa_rerun"
    if (liba / "RELEASE_VERIFICATION.json").is_file():
        target = records / "liba6656-gfa-rerun"
        target.mkdir(parents=True, exist_ok=True)
        verification = json.loads((liba / "RELEASE_VERIFICATION.json").read_text())
        shutil.copyfile(liba / "RELEASE_VERIFICATION.json", target / "RELEASE_VERIFICATION.json")
        first = sorted(liba.glob("run_*_threads"),
                       key=lambda path: int(path.name.split("_")[1]))[0] / "output"
        relative_report(first / "locus_recon_report_tcdB.tsv", first,
                        target / "standard_workflow_report_tcdB.tsv")
        sample = first / "LIBA6656_ST154"
        shutil.copyfile(sample / "LIBA6656_ST154_tcdB_qc_report.txt",
                        target / "standard_workflow_qc_report_tcdB.txt")
        shutil.copyfile(sample / "LIBA6656_ST154_tcdB_reconstructed.fasta",
                        target / "standard_workflow_reconstructed_tcdB.fasta")
        updated = {
            "standard_workflow_report": "standard_workflow_report_tcdB.tsv",
            "standard_workflow_qc_report": "standard_workflow_qc_report_tcdB.txt",
            "standard_workflow_reconstructed": "standard_workflow_reconstructed_tcdB.fasta",
        }
        checksum_rows = read_table(old / "liba6656-gfa-rerun" / "INPUT_CHECKSUMS.tsv")
        for row in checksum_rows:
            if row["role"] in updated:
                row["md5"] = digest_file(target / updated[row["role"]], "md5")
                row["sha256"] = digest_file(target / updated[row["role"]])
        write_table(target / "INPUT_CHECKSUMS.tsv", list(checksum_rows[0]),
                    [list(row.values()) for row in checksum_rows])
        readme = (old / "liba6656-gfa-rerun" / "README.md").read_text()
        section = args.liba_readme_section
        if section and section.is_file():
            start = readme.index("## Verification with the release code")
            end = readme.index("## Files")
            readme = readme[:start] + section.read_text().rstrip("\n") + "\n\n" + readme[end:]
            (target / "README.md").write_text(readme)
        print("liba6656-gfa-rerun: runs agree", verification["runs_agree"])
        built.append("liba6656-gfa-rerun")

    # frame and graph checks
    source = work / "frame_graph" / "frame-and-graph-evidence"
    if (source / "VALIDATION_SUMMARY.json").is_file():
        target = records / "frame-and-graph-evidence"
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(source, target)
        shutil.copyfile(old / "frame-and-graph-evidence" / "README.md", target / "README.md")
        built.append("frame-and-graph-evidence")

    # comparison with ARIBA and KMA: files the scoring does not write
    target = records / "comparator-benchmark"
    if (target / "comparator_benchmark_results.tsv").is_file():
        for name in ("raw_tool_outputs.tar.gz", "TOOL_VERSIONS.txt", "PROVENANCE.json",
                     "conda_environment.yml", "conda_list.txt", "kma_B7_9_gap_read_check.tsv",
                     "liba6656_analyst_built_candidates.fasta"):
            shutil.copyfile(old / "comparator-benchmark" / name, target / name)
        built.append("comparator-benchmark")

    for path in list(records.rglob("*")):
        if path.is_file() and (path.suffix in (".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq",
                                               ".ntf", ".nto") or path.name == ".DS_Store"):
            path.unlink()

    # README files whose text names the release or the source of a row: the
    # passages listed in readme_edits.py are replaced, each exactly once.
    try:
        from readme_edits import EDITS
    except ImportError:
        EDITS = {}
    for relative, replacements in EDITS.items():
        if not (records / relative).parent.is_dir():
            continue
        text = (old / relative).read_text()
        for before, after in replacements:
            if text.count(before) != 1:
                raise SystemExit(f"{relative}: passage found {text.count(before)} time(s), "
                                 f"expected once:\n{before}")
            text = text.replace(before, after)
        (records / relative).write_text(text)
        print(f"{relative}: {len(replacements)} passage(s) rewritten")
    print("records written:", ", ".join(built))


def legacy_fields_unchanged(earlier_table, rows, repo):
    """The comparison of the low-copy runner, without stopping at a difference."""
    sys.path.insert(0, str(repo))
    from validation.workflows.run_low_copy_23S import LEGACY_RESULT_FIELDS
    earlier = {(row["strain"], row["locus"]): row for row in read_table(earlier_table)}
    differences = []
    for row in rows:
        key = (row["strain"], row["locus"])
        for field in LEGACY_RESULT_FIELDS:
            if earlier.get(key, {}).get(field) != row.get(field):
                differences.append({
                    "strain": key[0], "locus": key[1], "field": field,
                    "release_1_0_0": earlier.get(key, {}).get(field),
                    "regenerated": row.get(field)})
    return {"unchanged": not differences, "differences": differences}


# ------------------------------------------------------------------ compare

VOLATILE = re.compile(
    r"(^|[._\[\]/ ])(run_date|elapsed_sec|seconds|seconds_min|seconds_max|created_utc|"
    r"generated_utc|completed_utc|started_at|tool_version|version|locus_recon_version|"
    r"source_sha256|git_commit|git_tree_clean|allele_file|message|python|biopython)"
    r"($|[._\[\]/ ])")
FLAG_COLUMNS = ("qc_flags", "dosage_flags", "consensus_flags", "graph_context_flags",
                "span_flags")
NEW_LABELS = ("MIXED_SITE_BELOW_MIXTURE_CALL", "DOSAGE_")


def flag_names(text):
    names, depth, current = [], 0, []
    for character in str(text):
        depth += {"(": 1, ")": -1}.get(character, 0)
        if character == ";" and depth == 0:
            names.append("".join(current))
            current = []
        else:
            current.append(character)
    names.append("".join(current))
    return ";".join(sorted(name.split("(")[0].strip() for name in names if name.strip()))


def flatten(value, prefix=""):
    if isinstance(value, dict):
        for key in sorted(value):
            yield from flatten(value[key], f"{prefix}.{key}" if prefix else str(key))
    elif isinstance(value, list):
        if all(not isinstance(item, (dict, list)) for item in value):
            yield prefix, json.dumps(value)
        else:
            for index, item in enumerate(value):
                yield from flatten(item, f"{prefix}[{index}]")
    else:
        yield prefix, value if isinstance(value, str) else json.dumps(value)


def keyed(rows, keys):
    """Rows by key columns; by row number when the key is absent or repeated."""
    if rows and keys and all(key in rows[0] for key in keys):
        names = [" | ".join(row[key] for key in keys) for row in rows]
        if len(set(names)) == len(names):
            return dict(zip(names, rows))
    return {f"row {index + 1}": row for index, row in enumerate(rows)}


class Comparison:
    def __init__(self):
        self.rows = []

    def add(self, analysis, record, case, quantity, old, new):
        old = "" if old is None else str(old)
        new = "" if new is None else str(new)
        kind = "run metadata" if VOLATILE.search(quantity) else "result"
        equal = "yes" if old == new else "no"
        difference = ""
        if equal == "no":
            try:
                difference = f"{abs(float(new) - float(old)):.6g}"
            except ValueError:
                pass
        self.rows.append([analysis, record, case, quantity, kind, old, new, equal, difference])

    def table(self, analysis, record, old_path, new_path, keys=(), rename=None):
        if not old_path.is_file() or not new_path.is_file():
            self.add(analysis, record, "", "file present",
                     "yes" if old_path.is_file() else "no",
                     "yes" if new_path.is_file() else "no")
            return
        old_rows, new_rows = read_table(old_path), read_table(new_path)
        if rename:
            for rows in (old_rows, new_rows):
                for row in rows:
                    rename(row)
        old, new = keyed(old_rows, keys), keyed(new_rows, keys)
        self.add(analysis, record, "", "number of rows", len(old_rows), len(new_rows))
        for case in list(old) + [case for case in new if case not in old]:
            first, second = old.get(case), new.get(case)
            if first is None or second is None:
                self.add(analysis, record, case, "row present",
                         "yes" if first else "no", "yes" if second else "no")
                continue
            for column in list(first) + [column for column in second if column not in first]:
                before, after = first.get(column), second.get(column)
                if column in keys:
                    continue
                if column not in first or column not in second:
                    self.add(analysis, record, case, f"{column} (column present)",
                             "yes" if column in first else "no",
                             "yes" if column in second else "no")
                    continue
                self.add(analysis, record, case, column, before, after)
                if column in FLAG_COLUMNS:
                    self.add(analysis, record, case, f"{column}: flag names",
                             flag_names(before), flag_names(after))
                if column == "message":
                    strip = lambda text: re.sub(r"\s*\|\s*[\d.]+s\s*$", "", text or "")
                    self.add(analysis, record, case, "summary line without the run time",
                             strip(before), strip(after))

    def json(self, analysis, record, old_path, new_path, subtree=None):
        if not old_path.is_file() or not new_path.is_file():
            self.add(analysis, record, "", "file present",
                     "yes" if old_path.is_file() else "no",
                     "yes" if new_path.is_file() else "no")
            return
        old, new = json.loads(old_path.read_text()), json.loads(new_path.read_text())
        if subtree:
            old, new = old.get(subtree, {}), new.get(subtree, {})
        first, second = dict(flatten(old)), dict(flatten(new))
        for key in list(first) + [key for key in second if key not in first]:
            if key not in first or key not in second:
                self.add(analysis, record, key, "key present",
                         "yes" if key in first else "no", "yes" if key in second else "no")
            else:
                self.add(analysis, record, "", key, first[key], second[key])

    def file(self, analysis, record, old_path, new_path):
        present = (old_path.is_file(), new_path.is_file())
        if not all(present):
            self.add(analysis, record, "", "file present",
                     "yes" if present[0] else "no", "yes" if present[1] else "no")
            return
        self.add(analysis, record, "", "file content (SHA-256 digest)",
                 digest_file(old_path), digest_file(new_path))

    def fasta(self, analysis, record, old_path, new_path, name=lambda text: text.split()[0]):
        if not old_path.is_file() or not new_path.is_file():
            self.file(analysis, record, old_path, new_path)
            return
        old = {name(header): sequence for header, sequence in read_fasta(old_path)}
        new = {name(header): sequence for header, sequence in read_fasta(new_path)}
        self.add(analysis, record, "", "number of sequences", len(old), len(new))
        for key in list(old) + [key for key in new if key not in old]:
            if key not in old or key not in new:
                self.add(analysis, record, key, "sequence present",
                         "yes" if key in old else "no", "yes" if key in new else "no")
                continue
            self.add(analysis, record, key, "sequence length", len(old[key]), len(new[key]))
            self.add(analysis, record, key, "sequence (SHA-256 digest)",
                     hashlib.sha256(old[key].encode()).hexdigest(),
                     hashlib.sha256(new[key].encode()).hexdigest())


def explain(row, notes):
    analysis, record, case, quantity, kind, old, new, equal, _ = row
    if equal == "yes":
        return ""
    if "copy_number_status" in quantity and new == "REVIEW":
        return ("label introduced in 1.0.1: copy_number_status is REVIEW when dosage_status "
                "is not ESTIMATED")
    if any(label in new and label not in old for label in NEW_LABELS):
        return "label introduced in 1.0.1"
    for pattern, text in notes:
        if pattern.search(" | ".join((analysis, record, case, quantity))):
            return text
    if kind == "run metadata":
        return ("run metadata (date, run time, version, digest, path or checksum of a "
                "regenerated file); expected to differ")
    return "UNEXPLAINED"


def command_compare(args):
    old, new = args.old, args.package / "records"
    comparison = Comparison()
    table, record_json, record_file, fasta = (
        comparison.table, comparison.json, comparison.file, comparison.fasta)

    def both(folder, name):
        return old / folder / name, new / folder / name

    name = "simulated benchmark"
    folder = "reference_results"
    table(name, f"{folder}/locus_recon_report_mockLocus.tsv",
          *both(folder, "locus_recon_report_mockLocus.tsv"), keys=("sample_id",))
    table(name, f"{folder}/validation_results.tsv", *both(folder, "validation_results.tsv"),
          keys=("sample_id",))
    record_json(name, f"{folder}/validation_summary.json",
                *both(folder, "validation_summary.json"))

    name, folder = "completeness benchmark", "completeness-benchmark"
    table(name, f"{folder}/completeness_benchmark_results.tsv",
          *both(folder, "completeness_benchmark_results.tsv"), keys=("case",))
    record_json(name, f"{folder}/completeness_benchmark_results.json",
                *both(folder, "completeness_benchmark_results.json"), subtree="summary")
    record_file(name, f"{folder}/benchmark_query.fasta", *both(folder, "benchmark_query.fasta"))
    for case in sorted((old / folder / "cases").glob("*")):
        for filename in ("assembly.fasta", "local_hits.tsv"):
            record_file(name, f"{folder}/cases/{case.name}/{filename}",
                        *both(folder, f"cases/{case.name}/{filename}"))

    name, folder = "tier calibration", "tier-calibration"
    for stage in ("benchmarks", "closed_genomes", "divergent_loci", "depth_series",
                  "depth_series_isolate"):
        table(name, f"{folder}/tier_calibration_{stage}.tsv",
              *both(folder, f"tier_calibration_{stage}.tsv"), keys=("case",))
    table(name, f"{folder}/tier_calibration_all_cases.tsv",
          *both(folder, "tier_calibration_all_cases.tsv"), keys=("stage", "case"))
    table(name, f"{folder}/locus_selection.tsv", *both(folder, "locus_selection.tsv"),
          keys=("gene",))
    record_json(name, f"{folder}/TIER_CALIBRATION.json", *both(folder, "TIER_CALIBRATION.json"))
    table(name, f"{folder}/INPUT_CHECKSUMS.tsv", *both(folder, "INPUT_CHECKSUMS.tsv"),
          keys=("path",))
    for filename in ("bait_23S_26695.fasta", "bait_cagA_26695.fasta", "bait_glmM_26695.fasta",
                     "bait_gyrB_26695.fasta", "bait_ureB_26695.fasta",
                     "closed_genome_manifest.tsv", "tier_calibration_tables.md"):
        record_file(name, f"{folder}/{filename}", *both(folder, filename))

    name, folder = "checks of single components", "frame-and-graph-evidence"
    record_json(name, f"{folder}/VALIDATION_SUMMARY.json",
                *both(folder, "VALIDATION_SUMMARY.json"))
    for filename in ("frame_orientation_invariance.tsv", "bait_validation_messages.tsv",
                     "length_profile_reporting.tsv", "placeholder_junction_cases.tsv",
                     "scoring_positive_control.tsv", "graph_path_scores.tsv"):
        table(name, f"{folder}/{filename}", *both(folder, filename))
    for filename in ("scoring_positive_control_candidates.fasta",
                     "scoring_positive_control_R1.fastq.gz",
                     "scoring_positive_control_R2.fastq.gz", "graph_path_scores_idxstats.txt"):
        record_file(name, f"{folder}/{filename}", *both(folder, filename))

    name, folder = "synthetic copy series", "constructed-copy-series"
    table(name, f"{folder}/results.tsv", *both(folder, "results.tsv"), keys=("case",))
    record_json(name, f"{folder}/evaluation.json", *both(folder, "evaluation.json"))
    record_json(name, f"{folder}/run_summary.json", *both(folder, "run_summary.json"))

    name, folder = "aphA1", "aphA1-copy-number"
    table(name, f"{folder}/results.tsv", *both(folder, "results.tsv"), keys=("run",))
    record_json(name, f"{folder}/evaluation.json", *both(folder, "evaluation.json"))
    record_json(name, f"{folder}/run_summary.json", *both(folder, "run_summary.json"))

    name, folder = "23S rRNA gene and gyrB", "low-copy-23S"
    table(name, f"{folder}/results.tsv", *both(folder, "results.tsv"), keys=("strain", "locus"))
    record_json(name, f"{folder}/evaluation.json", *both(folder, "evaluation.json"))
    record_json(name, f"{folder}/closed_truth_audit.json",
                *both(folder, "closed_truth_audit.json"))
    record_json(name, f"{folder}/run_summary.json", *both(folder, "run_summary.json"),
                subtree="panel")

    name, folder = "LIBA-6656 depth", "liba6656-depth"
    table(name, f"{folder}/LIBA6656_tcdB_depth_result.tsv",
          *both(folder, "LIBA6656_tcdB_depth_result.tsv"))
    table(name, f"{folder}/result_summary.tsv", *both(folder, "result_summary.tsv"))
    record_json(name, f"{folder}/evaluation.json", *both(folder, "evaluation.json"))
    folder = "depth-copy-number/overlap_null"
    table("LIBA-6656 overlap control", f"{folder}/null_dose_overlap_control.csv",
          *both(folder, "null_dose_overlap_control.csv"))
    record_json("LIBA-6656 overlap control", f"{folder}/null_dose_overlap_control.json",
                *both(folder, "null_dose_overlap_control.json"))

    name, folder = "LIBA-6656 standard workflow", "liba6656-gfa-rerun"
    record_json(name, f"{folder}/RELEASE_VERIFICATION.json",
                *both(folder, "RELEASE_VERIFICATION.json"))
    table(name, f"{folder}/standard_workflow_report_tcdB.tsv",
          *both(folder, "standard_workflow_report_tcdB.tsv"), keys=("sample_id",))
    fasta(name, f"{folder}/standard_workflow_reconstructed_tcdB.fasta",
          *both(folder, "standard_workflow_reconstructed_tcdB.fasta"))
    record_file(name, f"{folder}/standard_workflow_qc_report_tcdB.txt",
                *both(folder, "standard_workflow_qc_report_tcdB.txt"))
    table(name, f"{folder}/INPUT_CHECKSUMS.tsv", *both(folder, "INPUT_CHECKSUMS.tsv"),
          keys=("role",))

    name, folder = "comparison with ARIBA and KMA", "comparator-benchmark"

    def one_label(row):
        if row.get("tool") == "Locus-Recon":
            row["configuration_label"] = row["configuration"]
            row["configuration"] = "Locus-Recon"

    table(name, f"{folder}/comparator_benchmark_results.tsv",
          *both(folder, "comparator_benchmark_results.tsv"),
          keys=("dataset", "sample", "locus", "tool", "configuration"), rename=one_label)
    table(name, f"{folder}/comparator_benchmark_summary.tsv",
          *both(folder, "comparator_benchmark_summary.tsv"),
          keys=("dataset", "group", "tool", "configuration"), rename=one_label)
    fasta(name, f"{folder}/comparator_returned_sequences.fasta",
          *both(folder, "comparator_returned_sequences.fasta"))
    for filename in ("comparator_truth_sequences.fasta", "hpylori_catalogue.fasta",
                     "hpylori_catalogue_ariba_metadata.tsv", "liba6656_catalogue.fasta",
                     "liba6656_catalogue_ariba_metadata.tsv", "mock_catalogue.fasta",
                     "mock_catalogue_ariba_metadata.tsv"):
        record_file(name, f"{folder}/{filename}", *both(folder, filename))

    notes = []
    if args.notes and args.notes.is_file():
        for row in read_table(args.notes):
            notes.append((re.compile(row["pattern"]), row["explanation"]))
    header = ["analysis", "record", "case", "quantity", "kind", "value_1.0.0",
              "value_1.0.1", "equal", "absolute_difference", "explanation"]
    rows = [row + [explain(row, notes)] for row in comparison.rows]
    write_table(args.package / "old_vs_new.tsv", header, rows)

    counts = {}
    for row in rows:
        entry = counts.setdefault(row[0], {"values": 0, "equal": 0, "metadata": 0,
                                           "labels": 0, "explained": 0, "unexplained": 0})
        entry["values"] += 1
        if row[7] == "yes":
            entry["equal"] += 1
        elif row[9] == "UNEXPLAINED":
            entry["unexplained"] += 1
        elif row[9].startswith("run metadata"):
            entry["metadata"] += 1
        elif row[9].startswith("label introduced"):
            entry["labels"] += 1
        else:
            entry["explained"] += 1
    write_table(args.package / "old_vs_new_summary.tsv", [
        "analysis", "values_compared", "equal", "different_run_metadata",
        "different_labels_introduced_in_1.0.1", "different_explained_otherwise",
        "different_unexplained"],
        [[name, *entry.values()] for name, entry in counts.items()])
    for name, entry in counts.items():
        print(name, entry)
    print("UNEXPLAINED (first 60):")
    for row in [row for row in rows if row[9] == "UNEXPLAINED"][:60]:
        print("  ", " | ".join(str(value)[:110] for value in row[:7]))


# ------------------------------------------------------------------ summary

def command_summary(args):
    """Write summary.tsv: the headline numbers, each with the file it is read from."""
    package = args.package
    records, evidence = package / "records", package / "evidence"
    out = []

    def add(topic, quantity, value, source):
        out.append([topic, quantity, value, source])

    add("software", "Locus-Recon version", args.version, "PROVENANCE.json")
    add("software", "package source digest (SHA-256)", args.digest, "PROVENANCE.json")
    add("software", "git commit of the release tree", args.commit, "PROVENANCE.json")

    source = "old_vs_new_summary.tsv"
    rows = read_table(package / source)
    totals = {key: sum(int(row[key]) for row in rows) for key in rows[0] if key != "analysis"}
    add("comparison with the 1.0.0 records", "analyses compared", len(rows), source)
    for key, value in totals.items():
        add("comparison with the 1.0.0 records", key.replace("_", " "), value, source)
    for row in rows:
        different = int(row["values_compared"]) - int(row["equal"])
        add("comparison with the 1.0.0 records", f"{row['analysis']}: values compared / different",
            f"{row['values_compared']} / {different}", source)
    source = "old_vs_new.tsv"
    compared = read_table(package / source)
    results = [row for row in compared if row["kind"] == "result" and row["equal"] == "no"]
    add("comparison with the 1.0.0 records", "result values that differ (not run metadata)",
        len(results), source)
    add("comparison with the 1.0.0 records",
        "values carrying a label introduced in 1.0.1",
        sum(row["explanation"].startswith("label introduced") for row in compared), source)
    add("comparison with the 1.0.0 records",
        "values with copy_number_status REVIEW or a DOSAGE_ flag in the regenerated records",
        sum(("copy_number_status" in row["quantity"] and row["value_1.0.1"] == "REVIEW")
            or "DOSAGE_" in row["value_1.0.1"] for row in compared), source)
    flagged = sorted({row["case"] for row in compared
                      if "MIXED_SITE_BELOW_MIXTURE_CALL" in row["value_1.0.1"]
                      and row["record"].endswith("tier_calibration_closed_genomes.tsv")})
    add("comparison with the 1.0.0 records", "cases with MIXED_SITE_BELOW_MIXTURE_CALL",
        f"{len(flagged)} ({'; '.join(flagged)})", source)

    source = "records/reference_results/validation_summary.json"
    summary = json.loads((package / source).read_text())
    add("simulated benchmark", "core criteria met",
        f"{sum(summary['criteria'].values())} of {len(summary['criteria'])}", source)
    source = "records/reference_results/validation_results.tsv"
    rows = read_table(package / source)
    add("simulated benchmark", "samples", len(rows), source)
    add("simulated benchmark", "samples with a candidate identical to the truth",
        sum(row["exact_truth_match"] == "True" for row in rows), source)
    source = "records/reference_results/locus_recon_report_mockLocus.tsv"
    rows = read_table(package / source)
    add("simulated benchmark", "samples accepted (PASS)",
        sum(row["result_disposition"] == "PASS" for row in rows), source)
    source = "evidence/mock_run_and_repeat_checksums.tsv"
    rows = read_table(package / source)
    add("simulated benchmark", "kept files with the same bytes in a second run",
        f"{sum(row['identical'] == 'yes' for row in rows)} of {len(rows)}", source)

    source = "records/completeness-benchmark/completeness_benchmark_results.json"
    summary = json.loads((package / source).read_text())["summary"]
    add("completeness benchmark", "cases with the loss measured exactly",
        f"{summary['cases_exact']} of {summary['cases']}", source)

    source = "records/tier-calibration/TIER_CALIBRATION.json"
    overall = json.loads((package / source).read_text())["combined"]["overall"]
    add("tier calibration", "cases", overall["cases"], source)
    add("tier calibration", "cases accepted (tier HIGH)", overall["reached_top_tier"], source)
    add("tier calibration", "accepted cases identical to the truth", overall["top_tier_exact"],
        source)
    add("tier calibration", "accepted cases not identical to the truth",
        overall["false_accepts"], source)
    add("tier calibration", "exact 95% interval of the precision of accepted cases",
        "-".join(str(value) for value in overall["top_tier_precision_ci95"]), source)
    add("tier calibration", "cases not accepted but identical to the truth",
        overall["withheld_but_exact"], source)

    source = "exactness_check.tsv"
    rows = read_table(package / source)
    add("check of exactness against sequences", "cases with a kept candidate or a recorded absence",
        len(rows), source)
    add("check of exactness against sequences", "candidates identical to the truth sequence",
        sum(row["identical_to_truth"] == "yes" for row in rows), source)
    add("check of exactness against sequences", "cases where the check disagrees with the label of the runner",
        sum(row["check_agrees_with_label"] == "no" for row in rows), source)

    source = "records/comparator-benchmark/comparator_benchmark_results.tsv"
    rows = read_table(package / source)
    for tool, configuration in (("Locus-Recon", None), ("ARIBA", "default"), ("KMA", "default")):
        subset = [row for row in rows if row["dataset"] == "hpylori" and row["tool"] == tool
                  and (configuration is None or row["configuration"] == configuration)]
        add("comparison with ARIBA and KMA",
            f"H. pylori loci with a sequence identical to the truth, {tool}",
            f"{sum(row['comparison'] == 'identical' for row in subset)} of {len(subset)}", source)
    source = "old_vs_new.tsv"
    comparator = [row for row in compared if row["analysis"].startswith("comparison with ARIBA")
                  and row["record"].endswith("comparator_benchmark_results.tsv")]
    others = [row for row in comparator if " | Locus-Recon | " not in row["case"] + " | "]
    add("comparison with ARIBA and KMA", "values of ARIBA and KMA rows compared / different",
        f"{len(others)} / {sum(row['equal'] == 'no' for row in others)}", source)

    source = "records/constructed-copy-series/evaluation.json"
    evaluation = json.loads((package / source).read_text())
    add("synthetic copy series", "data sets", evaluation["cases"], source)
    add("synthetic copy series", "median absolute error of the dosage (copies)",
        evaluation["median_absolute_error"], source)
    add("synthetic copy series", "one-copy data sets called multicopy",
        evaluation["one_copy_false_multicopy"], source)
    add("synthetic copy series", "two- and three-copy data sets not called multicopy",
        evaluation["two_or_three_copy_missed"], source)

    source = "records/aphA1-copy-number/results.tsv"
    for row in read_table(package / source):
        add("aphA1", f"{row['strain']} ({row['run']}): dosage (ratio interval), call",
            f"{row['dosage_estimate']} ({row['ratio_ci_low']}-{row['ratio_ci_high']}), "
            f"{row['depth_call']}", source)

    source = "records/low-copy-23S/results.tsv"
    rows = read_table(package / source)
    for row in rows:
        add("23S rRNA gene and gyrB", f"{row['strain']} {row['locus']}: ratio (interval), call, "
            "graph count, consensus status",
            f"{row['ratio']} ({row['ratio_ci_low']}-{row['ratio_ci_high']}), {row['depth_call']}, "
            f"{row['graph_context_count'] or 'lower bound ' + row['graph_context_lower_bound']}, "
            f"{row['copy_number_status']}", source)
    review = [f"{row['strain']} {row['locus']}" for row in rows
              if row["dosage_status"] != "ESTIMATED"]
    supported = [f"{row['strain']} {row['locus']}" for row in rows
                 if row["dosage_status"] != "ESTIMATED" and row["copy_number_status"] == "SUPPORTED"]
    add("23S rRNA gene and gyrB", "loci whose dosage_status is not ESTIMATED",
        f"{len(review)} ({'; '.join(review) or 'none'})", source)
    add("23S rRNA gene and gyrB",
        "of these, loci whose copy_number_status is SUPPORTED",
        f"{len(supported)} ({'; '.join(supported) or 'none'})", source)

    source = "records/liba6656-depth/result_summary.tsv"
    row = read_table(package / source)[0]
    add("LIBA-6656 tcdB", "depth ratio (interval), call",
        f"{row['ratio']} ({row['ratio_ci_low']}-{row['ratio_ci_high']}), {row['depth_call']}",
        source)
    add("LIBA-6656 tcdB", "query-weighted dosage (interval), status",
        f"{row['dosage_estimate']} ({row['dosage_ci_low']}-{row['dosage_ci_high']}), "
        f"{row['dosage_status']}", source)
    add("LIBA-6656 tcdB", "dosage expected from the geometry alone",
        row["geometry_expected_dosage"], source)
    source = "records/depth-copy-number/overlap_null/null_dose_overlap_control.json"
    if (package / source).is_file():
        null = json.loads((package / source).read_text())
        add("LIBA-6656 tcdB", "control windows with two or more intervals",
            f"{null['n_multi_region']} of {null['n_windows_tried']}", source)
        add("LIBA-6656 tcdB", "dosage predicted by the control at the overlap of tcdB (95% interval)",
            f"{null['null_prediction_at_tcdB_overlap']} "
            f"({null['null_prediction_interval_95'][0]}-{null['null_prediction_interval_95'][1]})",
            source)
    source = "records/liba6656-gfa-rerun/RELEASE_VERIFICATION.json"
    verification = json.loads((package / source).read_text())
    report = verification["runs"][0]["report"]
    add("LIBA-6656 tcdB", "standard workflow: length, tier, disposition, mixed sites",
        f"{report['allele_length']} bp, {report['sequence_confidence']}, "
        f"{report['result_disposition']}, {report['mixed_site_count']}", source)
    add("LIBA-6656 tcdB", "threads of the two runs; the runs agree",
        f"{' and '.join(str(run['parameters']['threads']) for run in verification['runs'])}; "
        f"{str(verification['runs_agree']).lower()}", source)
    add("LIBA-6656 tcdB", "graph equal to the deposited graph (segments, links)",
        f"{str(verification['graph']['same_segment_sequences']).lower()}, "
        f"{str(verification['graph']['same_links_with_segments_matched_by_sequence']).lower()}",
        source)
    add("LIBA-6656 tcdB", "terminal paths equal to the deposited paths",
        f"{str(verification['terminal_paths']['same_sequence_under_same_name']).lower()} "
        f"({verification['terminal_paths']['paths']['rerun']} paths)", source)

    source = "evidence/frame_graph/thread_count_comparison.tsv"
    if (package / source).is_file():
        for row in read_table(package / source):
            if row["file"] == "graph_path_scores.tsv":
                add("checks of single components",
                    f"path scores, {row['run_a']} against {row['run_b']}: rows that differ "
                    "(largest difference in mapped reads)",
                    f"{row['rows_differing']} of {row['rows_compared']} "
                    f"({row['largest_difference_in_mapped_reads']})", source)

    source = "evidence/hpylori/draft_assemblies.tsv"
    if (package / source).is_file():
        rows = [row for row in read_table(package / source)
                if row["earlier_contigs_fasta_sha256"]]
        add("SPAdes drafts", "drafts with the same bytes as the drafts of September 2026 on the same host",
            f"{sum(row['byte_identical_to_earlier'] == 'yes' for row in rows)} of {len(rows)}",
            source)

    source = "block_size_sensitivity.tsv"
    rows = read_table(package / source)
    for dataset in sorted({row["dataset"] for row in rows}):
        for block in ("100", "250", "500", "1000"):
            subset = [row for row in rows if row["dataset"] == dataset and row["block_bp"] == block]
            widths = [float(row["ratio_ci_width"]) for row in subset]
            judged = [row for row in subset if row["ratio_ci_contains_reference"]]
            add("block size of the bootstrap",
                f"{dataset}, blocks of {block} bp: locus blocks; median interval width; "
                "intervals that contain the reference value",
                f"{'/'.join(sorted({row['bootstrap_locus_blocks'] for row in subset}, key=int))}; "
                f"{statistics.median(widths):.3f}; "
                f"{sum(row['ratio_ci_contains_reference'] == 'yes' for row in judged)} of "
                f"{len(judged)}", source)
        ratios = []
        for case in sorted({row["case"] for row in rows if row["dataset"] == dataset}):
            widths = [float(row["ratio_ci_width"]) for row in rows
                      if row["dataset"] == dataset and row["case"] == case]
            ratios.append(max(widths) / min(widths))
        add("block size of the bootstrap",
            f"{dataset}: widest over narrowest interval across the four block sizes, "
            "range over cases", f"{min(ratios):.2f}-{max(ratios):.2f}", source)

    for path in sorted((package / "tests").glob("*.pytest.txt")):
        last = path.read_text(errors="replace").strip().splitlines()[-1].strip("= ")
        add("unit tests", path.name.replace(".pytest.txt", "").replace("_", " "), last,
            f"tests/{path.name}")
    write_table(package / "summary.tsv", ["topic", "quantity", "value", "source_file"], out)
    print(f"summary.tsv: {len(out)} rows")


# -------------------------------------------------------------- integration

def command_integration(args):
    """Write INTEGRATION_RUN.md from the files of the simulated-benchmark run."""
    package, regen = args.package, args.regen
    evidence = package / "evidence"
    report = read_table(evidence / "mock" / "locus_recon_report_mockLocus.tsv")
    checked = {row["sample_id"]: row
               for row in read_table(evidence / "mock" / "validation" / "validation_results.tsv")}
    summary = json.loads(
        (evidence / "mock" / "validation" / "validation_summary.json").read_text())
    checksums = read_table(evidence / "mock_run_and_repeat_checksums.tsv")
    by_file = {row["file"]: row for row in checksums}
    tools = dict(line.split("\t", 1)
                 for line in (regen / "TOOL_VERSIONS.tsv").read_text().splitlines())
    steps = {}
    for line in (regen / "logs" / "steps.tsv").read_text().splitlines():
        fields = line.split("\t")
        steps[fields[0]] = dict(field.split("=", 1) for field in fields[1:])
    command = ""
    for line in (regen / "logs" / "mock.log").read_text(errors="replace").splitlines():
        if line.startswith("Running: "):
            command = line[len("Running: "):]
    command = command.replace(str(args.repo.resolve()) + "/", "").replace(
        str((regen / "work" / "mock").resolve()) + "/", "work/")
    command = re.sub(r"^\S*/python\b", "python", command)
    tests, total_tests = "", "unit"
    for path in sorted((package / "tests").glob("*.pytest.txt")):
        last = path.read_text(errors="replace").strip().splitlines()[-1]
        tests += f"| `{path.name}` | {last.strip('= ')} |\n"
        counts = [int(value) for value in re.findall(r"(\d+) (?:passed|failed|error)", last)]
        if counts:
            total_tests = str(sum(counts))

    def kind(name):
        base = name.split("/")[-1]
        return base.split("_mockLocus_")[-1] if "_mockLocus_" in base else base

    kinds = {}
    for row in checksums:
        entry = kinds.setdefault(kind(row["file"]), [0, 0])
        entry[0 if row["identical"] == "yes" else 1] += 1
    identical = sum(row["identical"] == "yes" for row in checksums)
    differing_samples = sorted({
        row["file"].split("/")[0] for row in checksums
        if row["identical"] == "no" and row["file"].endswith(
            ("contigs.fasta", "scaffolds.fasta", "assembly_graph_after_simplification.gfa",
             "allele_region.bed", "blast_step5.tsv"))})

    lines = [
        "# Integration run with the external tools", "",
        f"The {total_tests} unit tests under `tests/` replace SPAdes, BWA-MEM2, samtools and "
        "BLAST+ with",
        "stubs. They check the logic of the package and never start those programs. This file",
        "documents one run of the complete workflow with the real programs, on the smallest",
        "data set of the repository: the 13 simulated samples of `validation/mock_dataset/`.",
        "The run can be repeated in about two minutes, and its outputs compared with the",
        "checksums below.", "",
        "## Environment of the recorded run", "",
        "| Component | Version |", "|---|---|",
    ]
    for name in ("locus-recon", "python", "spades", "bwa-mem2", "samtools", "blastn",
                 "makeblastdb", "biopython", "numpy", "platform", "cpu"):
        if name in tools:
            lines.append(f"| {name} | {tools[name]} |")
    lines += [
        "", f"Package source digest (`source_sha256` of `validation/common.py`): `{args.digest}`.",
        f"Git commit of the tree: `{args.commit}`.", "",
        "## Commands", "",
        "```bash",
        "mamba env create -f environment.yml",
        "mamba activate locus-recon",
        "python -m pip install -e .",
        "python validation/run_mock_validation.py \\",
        "    --results-dir work/mock_results --threads 2 --memory-gb 4 --force",
        "```", "",
        "The runner checks the data set, starts the workflow with the command below and then",
        "evaluates the result against the simulated truth "
        "(`validation/evaluate_mock_results.py`):", "",
        "```bash", command, "```", "",
        "For every sample the workflow runs `makeblastdb` and `blastn` (discovery), "
        "`bwa-mem2 index` and",
        "`bwa-mem2 mem` with `samtools` (mapping to the draft and recruitment of reads), "
        "`spades.py`",
        "(local assembly), `blastn` (location and validation of the candidate) and "
        "`bwa-mem2` with",
        "`samtools mpileup` (per-base support of the candidate).", "",
        "## Expected result", "",
        f"`validation_summary.json`: `all_core_criteria_passed` is "
        f"{str(summary['all_core_criteria_passed']).lower()}; "
        f"{sum(summary['criteria'].values())} of {len(summary['criteria'])} criteria are met.",
        "",
        "| Sample | Status | Tier | Disposition | Length (bp) | Identical to the truth | "
        "Mixed sites | SHA-256 of `<sample>_mockLocus_reconstructed.fasta` |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in report:
        sample = row["sample_id"]
        digest = by_file.get(f"{sample}/{sample}_mockLocus_reconstructed.fasta", {})
        lines.append(
            f"| `{sample}` | {row['workflow_status']} | {row['sequence_confidence'] or '—'} | "
            f"{row['result_disposition'] or '—'} | {row['allele_length'] or '—'} | "
            f"{checked[sample]['exact_truth_match'].lower()} | {row['mixed_site_count'] or '—'} | "
            f"{'`' + digest['sha256_run'] + '`' if digest else 'no sequence'} |")
    lines += [
        "", "Files of the results directory that summarise the run:", "",
        "| File | Bytes | SHA-256 | Same bytes in a second run |", "|---|---|---|---|",
    ]
    for name in ("mockLocus_accepted_alleles.fasta", "mockLocus_review_required_candidates.fasta",
                 "mockLocus_hold_candidates.fasta", "mockLocus_reconstructed_candidates.fasta",
                 "validation/validation_summary.json", "validation/validation_results.tsv",
                 "locus_recon_report_mockLocus.tsv"):
        if name in by_file:
            row = by_file[name]
            lines.append(f"| `{name}` | {row['size_bytes']} | `{row['sha256_run']}` | "
                         f"{row['identical']} |")
    lines += [
        "", "## What is reproducible to the byte", "",
        f"The command was run twice on the same machine "
        f"({steps.get('mock', {}).get('seconds', '?')} s and "
        f"{steps.get('mock_repeat', {}).get('seconds', '?')} s). Of the {len(checksums)} files "
        f"kept from the first run, {identical} had the same bytes in the second.", "",
        "| Kind of file | Same bytes | Different bytes |", "|---|---|---|",
    ]
    for name, (same, different) in sorted(kinds.items()):
        lines.append(f"| `{name}` | {same} | {different} |")
    lines += [
        "",
        "The candidate sequences, the QC reports, the per-base tables "
        "(`allele_remap.support.tsv`) and",
        "`validation_summary.json` are the files to compare. Logs, the batch report, "
        "`validation_results.tsv`,",
        "`run_manifest.json` and `workflow_provenance.json` carry dates, run times or absolute "
        "paths and differ",
        "between runs. SPAdes writes its own paths and times in `params.txt` and `spades.log`.",
    ]
    if differing_samples:
        lines += [
            f"In {len(differing_samples)} sample(s) ({', '.join('`' + name + '`' for name in differing_samples)}) "
            "the local SPAdes assembly itself differed between the two runs, and with it the "
            "files that describe the contig; the candidate sequence, its QC report and its "
            "per-base table were the same.",
        ]
    lines += [
        "", "To compare a new run with the recorded one:", "",
        "```bash",
        "cd work/mock_results",
        "sha256sum */*_mockLocus_reconstructed.fasta mockLocus_*.fasta "
        "validation/validation_summary.json",
        "```", "",
        "The complete list of checksums of the recorded run is "
        "`evidence/mock_run_and_repeat_checksums.tsv` in the",
        "regeneration package (`validation/regenerated-1.0.1/`), and the files themselves are "
        "under `evidence/mock/`.", "",
        "## Relation to the unit tests", "",
        "| Test run | Last line of pytest |", "|---|---|", tests.rstrip("\n"), "",
        "All external programs are replaced with stubs in those tests. A change that breaks "
        "the call of an",
        "external program, or that depends on its real output, is therefore visible only in "
        "a run such as this",
        "one. The other analyses under `validation/` are also runs with the real programs; "
        "they need public",
        "reads of 0.2 to 1 GB per sample and take hours.", "",
    ]
    (package / "records").mkdir(parents=True, exist_ok=True)
    (package / "records" / "INTEGRATION_RUN.md").write_text("\n".join(lines))
    print("INTEGRATION_RUN.md written;", len(lines), "lines")


# ----------------------------------------------------------------- manifest

def command_manifest(args):
    package, regen = args.package, args.regen
    tools = {}
    for line in (regen / "TOOL_VERSIONS.tsv").read_text().splitlines():
        name, _, version = line.partition("\t")
        tools[name] = version
    steps, commands = [], {}
    for line in (regen / "logs" / "steps.tsv").read_text().splitlines():
        fields = line.split("\t")
        steps.append({"step": fields[0], **dict(field.split("=", 1) for field in fields[1:])})
    for log in sorted((regen / "logs").glob("*.log")):
        found = [line.lstrip("+ ").rstrip() for line in log.read_text(errors="replace").splitlines()
                 if re.match(r"^\++ (python|spades\.py|locus-recon|tar|pip|conda activate)", line)]
        if found:
            commands[log.stem] = found
    inputs = []
    table = regen / "inputs" / "INPUT_SHA256.tsv"
    if table.is_file():
        for line in table.read_text().splitlines():
            path, size, digest = line.split("\t")
            inputs.append({"path": path, "size_bytes": int(size), "sha256": digest,
                           "source": "ENA run" if path.endswith(".fastq.gz")
                           else "draft assembly of LIBA-6656 supplied by the author"})
    for name in ("aphA1/results/workflow_provenance.json",
                 "low_copy_23S/results/workflow_provenance.json",
                 "liba6656_depth/results/workflow_provenance.json",
                 "constructed_copy_series/results/workflow_provenance.json"):
        path = args.work / name
        if path.is_file():
            record = json.loads(path.read_text())
            for entry in record["inputs"]:
                inputs.append({**entry, "source": f"inputs recorded by {record['workflow']}"})
    seeds = {}
    for path in sorted((package / "records").rglob("*.json")) + sorted(
            (args.repo / "validation" / "mock_dataset").glob("dataset_manifest.json")):
        try:
            payload = json.loads(path.read_text())
        except ValueError:
            continue
        for key, value in flatten(payload):
            if re.search(r"(^|[._])(seed|base_seed)$", key) and len(seeds) < 200:
                label = (path.relative_to(package).as_posix() if package in path.parents
                         else f"validation/mock_dataset/{path.name}")
                seeds[f"{label}:{key}"] = value
    write_json(package / "PROVENANCE.json", {
        "schema_version": "1.0",
        "workflow": "regeneration of the validation analyses",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "software": {"name": "locus-recon", "version": args.version,
                     "source_sha256": args.digest, "git_commit": args.commit},
        "tools": tools, "steps": steps, "commands": commands, "seeds": seeds,
        "inputs": inputs,
        "driver": "validation/regenerate/regenerate_validation.sh",
    })
    lines = []
    for path in sorted(package.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            lines.append(f"{digest_file(path)}  {path.relative_to(package).as_posix()}")
    (package / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")
    total = sum(path.stat().st_size for path in package.rglob("*") if path.is_file())
    print(f"manifest: {len(lines)} files, {total} bytes")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("evidence", "records", "compare", "integration",
                                            "summary", "manifest"))
    parser.add_argument("--regen", type=Path, required=True, help="working directory")
    parser.add_argument("--package", type=Path, required=True, help="package directory")
    parser.add_argument("--old", type=Path, help="validation/ of the earlier release")
    parser.add_argument("--repo", type=Path, help="working clone of the release")
    parser.add_argument("--earlier-tier-work", type=Path, default=None,
                        help="work directory of an earlier tier calibration on this host")
    parser.add_argument("--liba-readme-section", type=Path, default=None)
    parser.add_argument("--notes", type=Path, default=None,
                        help="TSV with columns pattern and explanation")
    parser.add_argument("--version", default="1.0.1")
    parser.add_argument("--digest", default="")
    parser.add_argument("--commit", default="")
    args = parser.parse_args()
    args.work = args.regen / "work"
    args.driver_records = args.regen / "records"
    args.package.mkdir(parents=True, exist_ok=True)
    {"evidence": command_evidence, "records": command_records,
     "compare": command_compare, "integration": command_integration,
     "summary": command_summary, "manifest": command_manifest}[args.command](args)


if __name__ == "__main__":
    main()
