#!/usr/bin/env python3
"""Evidence, records and comparison of the analyses rerun with Locus-Recon 1.0.2.

The functions that copy evidence files, read tables and compare two records
value by value are those of ``collect_regenerated.py``, the helper of the 1.0.1
regeneration, which must lie beside this file.  The records of 1.0.1 are the
files under ``validation/`` of the release tree; the evidence of 1.0.1 is read
from the package of that regeneration.

Commands, in the order of the driver: evidence, records, integration, compare,
finish, summary, manifest.
"""
from __future__ import annotations

import argparse
import difflib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import collect_regenerated as cr  # noqa: E402

RULE_COLUMNS = ("copy_number_status", "copy_number_lower_bound", "copy_number_upper_bound")
RULE_METHOD = "GRAPH_COUNT_OVER_DISCORDANT_DEPTH"
RULE_TEXT = (
    "rule changed in 1.0.2: an integer graph count of two or more that the depth-ratio "
    "interval excludes is returned with copy_number_status REVIEW, copy_number_lower_bound "
    "equal to the count and no upper bound (1.0.1: SUPPORTED, both bounds equal to the count)")
METADATA_TEXT = ("run metadata (date, run time, version, digest, path or checksum of a "
                 "regenerated file); expected to differ")
MOCK_KINDS = ("_reconstructed.fasta", "_qc_report.txt", "allele_remap.support.tsv",
              "_alleles.fasta", "_candidates.fasta")


def checksum_table(new_dir, other_dir, destination, labels):
    """Every file under new_dir against the file of the same name under other_dir."""
    rows = []
    for path in sorted(Path(new_dir).rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(new_dir)
        other = Path(other_dir) / relative
        first = cr.digest_file(path)
        second = cr.digest_file(other) if other.is_file() else ""
        rows.append([relative.as_posix(), path.stat().st_size, first, second,
                     "yes" if first == second else "no"])
    cr.write_table(destination, ["file", "size_bytes", labels[0], labels[1], "identical"], rows)
    return rows


# ------------------------------------------------------------------ evidence

def command_evidence(args):
    work, out = args.work, args.package / "evidence"
    if out.exists():
        shutil.rmtree(out)
    skipped = []
    previous = args.prev_package / "evidence"

    # simulated benchmark: the run that is the record, and its repetition
    mock = work / "mock" / "mock_results"
    cr.copy_run(mock, out / "mock", skipped)
    for name in ("validation_results.tsv", "validation_summary.json", "workflow_provenance.json"):
        cr.copy_file(mock / "validation" / name, out / "mock" / "validation" / name, skipped)
    checksum_table(out / "mock", work / "mock" / "mock_results_repeat",
                   out / "mock_run_and_repeat_checksums.tsv", ("sha256_run", "sha256_repeat"))
    repeat = work / "mock" / "mock_results_repeat"
    for name in ("locus_recon_report_mockLocus.tsv", "run_manifest.json"):
        cr.copy_file(repeat / name, out / "mock_repeat" / name, skipped)
    for name in ("validation_results.tsv", "validation_summary.json"):
        cr.copy_file(repeat / "validation" / name, out / "mock_repeat" / "validation" / name,
                     skipped)
    checksum_table(out / "mock", previous / "mock", out / "mock_1.0.2_vs_1.0.1_checksums.tsv",
                   ("sha256_1.0.2", "sha256_1.0.1"))

    # LIBA-6656: the two runs of the standard workflow
    liba = work / "liba6656_gfa_rerun"
    for run in sorted(liba.glob("run_*_threads")):
        cr.copy_run(run / "output", out / "liba6656" / run.name, skipped)
        cr.copy_file(run / "samples.tsv", out / "liba6656" / run.name / "samples.tsv", skipped)
    for name in ("tcdB_terminal_paths_8000_9000.tsv", "RELEASE_VERIFICATION.json"):
        cr.copy_file(liba / name, out / "liba6656" / name, skipped)
    paths = liba / "tcdB_terminal_paths_8000_9000.fasta"
    if paths.is_file():
        (out / "liba6656").mkdir(parents=True, exist_ok=True)
        (out / "liba6656" / "tcdB_terminal_paths_8000_9000.fasta.sha256").write_text(
            f"{cr.digest_file(paths)}  tcdB_terminal_paths_8000_9000.fasta\n")
    checksum_table(out / "liba6656", previous / "liba6656",
                   out / "liba6656_1.0.2_vs_1.0.1_checksums.tsv",
                   ("sha256_1.0.2", "sha256_1.0.1"))

    # 23S rRNA gene and gyrB: results, provenance, logs (no alignments)
    def copy_small(directory, target, patterns):
        for pattern in patterns:
            for path in sorted(Path(directory).glob(pattern)):
                cr.copy_file(path, out / "copy_number" / target / path.name, skipped)

    low = work / "low_copy_23S"
    copy_small(low / "results", "low_copy_23S", ("*.tsv", "*.json"))
    copy_small(low / "fixed_baits", "low_copy_23S/fixed_baits", ("*.fasta",))
    for strain in cr.STRAINS:
        copy_small(low / "samples" / strain, f"low_copy_23S/samples/{strain}",
                   ("fixed_baits_vs_draft.tsv", "*.mapping.log", "spades.command.txt"))
    reused = low / "drafts_reused.tsv"
    if reused.is_file():
        rows = [line.split("\t") for line in reused.read_text().splitlines() if line]
        cr.write_table(out / "copy_number" / "low_copy_23S" / "drafts_reused.tsv",
                       ["strain", "file", "sha256_used_by_1.0.2", "sha256_used_by_1.0.1",
                        "identical"],
                       [row + ["yes" if row[2] == row[3] else "no"] for row in rows])
    checksum_table(out / "copy_number" / "low_copy_23S", previous / "copy_number" / "low_copy_23S",
                   out / "low_copy_23S_1.0.2_vs_1.0.1_checksums.tsv",
                   ("sha256_1.0.2", "sha256_1.0.1"))

    # locus-recon-copy-number without a graph on the LIBA-6656 alignment
    copy_small(work / "liba6656_consensus_nograph", "liba6656_consensus_without_graph",
               ("consensus_without_graph.json", "consensus_without_graph.tsv", "stdout.txt",
                "stderr.txt", "command.txt"))

    cr.write_table(args.package / "FILES_LEFT_ON_HOST.tsv", ["path", "size_bytes", "reason"],
                   [[path, size, reason] for path, size, reason in sorted(skipped)])
    files = [path for path in out.rglob("*") if path.is_file()]
    print(f"evidence: {len(files)} files, {sum(path.stat().st_size for path in files)} bytes; "
          f"{len(skipped)} file(s) above 5 MB left on the host")


# ------------------------------------------------------------------- records

def long_date(text):
    moment = datetime.fromisoformat(text.replace("Z", "+00:00"))
    return f"{moment.day} {moment.strftime('%B %Y')}"


def replace_once(text, old, new):
    """Replace one passage, whatever its line breaks; fail unless it occurs once."""
    pattern = r"\s+".join(re.escape(word) for word in old.split())
    found = re.findall(pattern, text)
    if len(found) != 1:
        raise ValueError(f"passage found {len(found)} times: {old[:60]}")
    return re.sub(pattern, lambda _: new, text)


def liba_readme(args, target, verification):
    """The README of the LIBA-6656 record with its verification section for this release.

    The section describes the graph, the paths and the reported sequence.  It is
    rewritten only when those are verified equal to the 1.0.1 record; otherwise
    the reasons are returned and the README is left out of the records.
    """
    old = args.old / "liba6656-gfa-rerun"
    liba = args.work / "liba6656_gfa_rerun"
    previous_liba = args.previous / "work" / "liba6656_gfa_rerun"
    reasons = []
    before = dict(cr.flatten(json.loads((old / "RELEASE_VERIFICATION.json").read_text())))
    after = dict(cr.flatten(verification))
    for key in sorted(set(before) | set(after)):
        if before.get(key) != after.get(key) and not cr.VOLATILE.search(key):
            reasons.append(f"RELEASE_VERIFICATION.json {key}: {before.get(key)} -> {after.get(key)}")
    for name in ("standard_workflow_reconstructed_tcdB.fasta", "standard_workflow_qc_report_tcdB.txt"):
        if (old / name).read_bytes() != (target / name).read_bytes():
            reasons.append(f"{name} differs in bytes from the 1.0.1 record")
    (row_before,), (row_after,) = (cr.read_table(old / "standard_workflow_report_tcdB.tsv"),
                                   cr.read_table(target / "standard_workflow_report_tcdB.tsv"))
    strip = lambda text: re.sub(r"\s*\|\s*[\d.]+s\s*$", "", text or "")
    for column in row_before:
        if column in ("tool_version", "run_date", "elapsed_sec"):
            continue
        first, second = row_before[column], row_after.get(column)
        if column == "message":
            first, second = strip(first), strip(second)
        if first != second:
            reasons.append(f"report row, {column}: {first} -> {second}")
    graph = "output/LIBA6656_ST154/spades_local/assembly_graph_after_simplification.gfa"
    threads = sorted(int(path.name.split("_")[1]) for path in liba.glob("run_*_threads"))
    for count in threads:
        if cr.digest_file(liba / f"run_{count}_threads" / graph) != cr.digest_file(
                previous_liba / f"run_{count}_threads" / graph):
            reasons.append(f"graph of the run with {count} threads differs from that of the 1.0.1 run")
    for name in ("tcdB_terminal_paths_8000_9000.fasta", "tcdB_terminal_paths_8000_9000.tsv"):
        if cr.digest_file(liba / name) != cr.digest_file(previous_liba / name):
            reasons.append(f"{name} differs in bytes from the file of the 1.0.1 run")
    if threads != [4, 6]:
        reasons.append(f"thread counts of the runs are {threads}, the text says four and six")
    if reasons:
        return reasons
    created = next((value for key, value in after.items()
                    if re.search(r"(created|generated|completed)_utc$", key)), None)
    if created is None:
        created = datetime.fromtimestamp(
            (liba / "RELEASE_VERIFICATION.json").stat().st_mtime, timezone.utc).isoformat()
    text = (old / "README.md").read_text()
    start = text.index("## Verification with the release code")
    end = text.index("## Files")
    section = text[start:end]
    section = replace_once(
        section,
        "On 18 November 2026 the two commands above were rerun with Locus-Recon 1.0.1 "
        "(source digest `2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0`)",
        f"On {long_date(created)} the two commands above were rerun with Locus-Recon "
        f"{args.version}\n(source digest `{args.digest}`)")
    section = replace_once(
        section,
        "The same verification was made with release 1.0.0 on 28 September 2026, with four and "
        "eight threads; its record is in the history of the repository (tag `v1.0.0`), and "
        "`../regenerated-1.0.1/old_vs_new.tsv` lists both records value by value.",
        "The same\nverification was made with release 1.0.0 on 28 September 2026, with four and\n"
        "eight threads, and with release 1.0.1 on 18 November 2026, with four and six\n"
        "threads; their records are in the history of the repository (tags `v1.0.0` and\n"
        "`v1.0.1`). `../regenerated-1.0.1/old_vs_new.tsv` and\n"
        f"`../regenerated-{args.version}/old_vs_new_{args.previous_version}_vs_{args.version}.tsv` "
        "list the\nrecords value by value.")
    section = replace_once(
        section,
        "- The reported sequence and the QC report have the same bytes as those that release "
        "1.0.0 wrote. The report row differs from that of 1.0.0 in its version, run date and "
        "run time only.",
        "- The reported sequence and the QC report have the same bytes as those that\n"
        "  releases 1.0.0 and 1.0.1 wrote. The report row differs from theirs in its\n"
        "  version, run date and run time only.")
    section = replace_once(
        section,
        "and the three `standard_workflow_*` files are the output of release 1.0.1.",
        f"and the three `standard_workflow_*` files are the output of release {args.version}.")
    (target / "README.md").write_text(text[:start] + section + text[end:])
    return []


def command_records(args):
    work, records, old = args.work, args.package / "records", args.old
    if records.exists():
        shutil.rmtree(records)
    records.mkdir(parents=True)

    # simulated benchmark: the three files written by the driver, and the text record
    shutil.copytree(args.regen / "records" / "reference_results", records / "reference_results")
    cr.mock_reproduction(records, work, args.package)

    # 23S rRNA gene and gyrB
    results = work / "low_copy_23S" / "results"
    target = records / "low-copy-23S"
    target.mkdir()
    shutil.copyfile(results / "low_copy_results.tsv", target / "results.tsv")
    for name in ("closed_truth_audit.json", "evaluation.json", "workflow_provenance.json"):
        if (results / name).is_file():
            shutil.copyfile(results / name, target / name)
    rows = cr.read_table(target / "results.tsv")
    rrna = [row for row in rows if row["locus"] == "23S"]
    controls = [row for row in rows if row["locus"] == "gyrB"]
    earlier = json.loads((old / "low-copy-23S" / "run_summary.json").read_text())
    baseline = work / "low_copy_23S" / "results_of_release_1.0.0.tsv"
    baseline.write_bytes(subprocess.run(
        ["git", "-C", str(args.repo), "show", "v1.0.0:validation/low-copy-23S/results.tsv"],
        check=True, capture_output=True).stdout)
    legacy = cr.legacy_fields_unchanged(baseline, rows, args.repo)
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
        by_rule = [row["strain"] for row in rrna
                   if row["copy_number_status"] == "REVIEW"
                   and row["copy_number_method"] == RULE_METHOD]
        interpretation = (
            f"The aggregate-depth rule classified {multicopy} of {len(rrna)} two-copy 23S "
            f"loci as MULTICOPY_DEPTH and {false_positive} of {len(controls)} single-copy "
            f"gyrB controls as MULTICOPY_DEPTH. Observed 23S depth ratios were "
            f"{min(ratios):.3f}-{max(ratios):.3f}. The copy-number consensus was an exact, "
            f"supported count of two for {panel['exact_two_copy_23S_consensus_calls']} of "
            f"{len(rrna)} 23S loci, a lower bound (copy_number_kind LOWER_BOUND) for "
            f"{panel['two_copy_23S_graph_lower_bounds']}, and an exact, supported count of "
            f"one for {panel['exact_single_copy_gyrB_consensus_calls']} of {len(controls)} "
            f"gyrB controls. 23S loci with copy_number_status REVIEW: "
            f"{', '.join(review) or 'none'}. ")
        if by_rule:
            interpretation += (
                f"In {len(by_rule)} of them ({', '.join(by_rule)}) the graph holds two flanking "
                "contexts and the depth-ratio interval excludes two; release "
                f"{args.version} returns such a count for review, with copy_number_lower_bound "
                "equal to the count and no upper bound, where release "
                f"{args.previous_version} returned SUPPORTED with both bounds equal to the "
                "count. ")
        interpretation += (
            f"all_primary_criteria_passed is "
            f"{str(panel['all_primary_criteria_passed']).lower()}. The panel informed the "
            "design of the graph method; it is a biological demonstration, not an "
            "independent performance estimate.")
        cr.write_json(target / "run_summary.json", {
            "schema_version": "1.0", "workflow": provenance["workflow"],
            "completed_utc": provenance["created_utc"],
            "software": {key: provenance["software"][key]
                         for key in ("name", "version", "source_sha256")},
            "panel": panel, "interpretation": interpretation,
            "legacy_depth_field_differences": legacy["differences"],
            "execution_note": (
                "The five SPAdes --isolate drafts are those of the regeneration with release "
                f"{args.previous_version} (assembled by the closed-genome stage of its tier "
                "calibration with 4 threads and a 14 GB memory cap, and reused by its 23S "
                "analysis); they were copied with unchanged checksums, and the workflow keeps "
                "an existing draft, so memory_gb below applied to no assembly in this run. "
                "The GenBank records of the closed genomes are the files that regeneration "
                "downloaded. The reads were mapped again and every later stage was run again. "
                "The GFA checksums are recorded in workflow_provenance.json."),
            "panel_of_release_1_0_1": earlier.get("panel"),
            "panel_of_release_1_0_0": earlier.get("panel_of_release_1_0_0"),
            "parameters": {key: value for key, value in provenance["parameters"].items()
                           if key != "selection"},
            "result_sha256": cr.digest_file(target / "results.tsv"),
            "evaluation_sha256": cr.digest_file(target / "evaluation.json"),
            "truth_audit_sha256": cr.digest_file(target / "closed_truth_audit.json"),
            "provenance_sha256": cr.digest_file(target / "workflow_provenance.json"),
        })
        print("low-copy-23S: legacy depth fields equal to release 1.0.0:", legacy["unchanged"],
              "| REVIEW:", review, "| by the rule of 1.0.2:", by_rule)
    else:
        print("low-copy-23S: INCOMPLETE, evaluation.json or workflow_provenance.json missing")

    # LIBA-6656, standard workflow and release verification
    liba = work / "liba6656_gfa_rerun"
    target = records / "liba6656-gfa-rerun"
    target.mkdir()
    verification = json.loads((liba / "RELEASE_VERIFICATION.json").read_text())
    shutil.copyfile(liba / "RELEASE_VERIFICATION.json", target / "RELEASE_VERIFICATION.json")
    first = sorted(liba.glob("run_*_threads"),
                   key=lambda path: int(path.name.split("_")[1]))[0] / "output"
    cr.relative_report(first / "locus_recon_report_tcdB.tsv", first,
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
    checksum_rows = cr.read_table(old / "liba6656-gfa-rerun" / "INPUT_CHECKSUMS.tsv")
    for row in checksum_rows:
        if row["role"] in updated:
            row["md5"] = cr.digest_file(target / updated[row["role"]], "md5")
            row["sha256"] = cr.digest_file(target / updated[row["role"]])
    cr.write_table(target / "INPUT_CHECKSUMS.tsv", list(checksum_rows[0]),
                   [list(row.values()) for row in checksum_rows])
    reasons = liba_readme(args, target, verification)
    print("liba6656-gfa-rerun: runs agree", verification["runs_agree"],
          "| software", verification.get("software"))
    if reasons:
        print("LIBA_README_NOT_UPDATED: the record differs from that of 1.0.1 in")
        for reason in reasons:
            print("   ", reason)
    else:
        print("LIBA_README_UPDATED: graph, paths, sequence, QC report and verification record "
              "equal those of the 1.0.1 run")

    # locus-recon-copy-number without a graph, LIBA-6656 tcdB
    target = records / "liba6656-depth"
    target.mkdir()
    for name in ("consensus_without_graph.json", "consensus_without_graph.tsv"):
        shutil.copyfile(work / "liba6656_consensus_nograph" / name, target / name)
    print("records:", sorted(path.relative_to(records).as_posix()
                             for path in records.rglob("*") if path.is_file()))


def command_integration(args):
    cr.command_integration(args)
    path = args.package / "records" / "INTEGRATION_RUN.md"
    text = path.read_text()
    assert text.count("validation/regenerated-1.0.1/") == 1
    path.write_text(text.replace("validation/regenerated-1.0.1/",
                                 f"validation/regenerated-{args.version}/"))
    print("INTEGRATION_RUN.md written;", len(text.splitlines()), "lines")


# ------------------------------------------------------------------- compare

def rule_cases(args):
    """Loci whose consensus the changed rule applies to, in either record."""
    cases = set()
    for table in (args.old / "low-copy-23S" / "results.tsv",
                  args.package / "records" / "low-copy-23S" / "results.tsv"):
        for row in cr.read_table(table):
            if row["copy_number_method"] == RULE_METHOD and row["depth_call"] == "MULTICOPY_DEPTH":
                cases.add((row["strain"], row["locus"]))
    return cases


def explain(row, notes, cases):
    analysis, record, case, quantity, kind, old, new, equal, _ = row
    if equal == "yes":
        return ""
    if record.endswith("low-copy-23S/results.tsv") and quantity in RULE_COLUMNS and any(
            strain in case and locus in case for strain, locus in cases):
        # only the change the release describes counts as the rule: the status
        # goes from SUPPORTED to REVIEW and the upper bound becomes empty; the
        # lower bound is expected to stay, so a change of it is not explained here
        if quantity == "copy_number_status" and (old, new) == ("SUPPORTED", "REVIEW"):
            return RULE_TEXT
        if quantity == "copy_number_upper_bound" and old != "" and new == "":
            return RULE_TEXT
    for pattern, text in notes:
        if pattern.search(" | ".join((analysis, record, case, quantity))):
            return text
    if kind == "run metadata":
        return METADATA_TEXT
    return "UNEXPLAINED"


def command_compare(args):
    old, new = args.old, args.package / "records"
    comparison = cr.Comparison()
    table, record_json, record_file, fasta = (
        comparison.table, comparison.json, comparison.file, comparison.fasta)

    def both(folder, name):
        return old / folder / name, new / folder / name

    name, folder = "simulated benchmark", "reference_results"
    table(name, f"{folder}/locus_recon_report_mockLocus.tsv",
          *both(folder, "locus_recon_report_mockLocus.tsv"), keys=("sample_id",))
    table(name, f"{folder}/validation_results.tsv", *both(folder, "validation_results.tsv"),
          keys=("sample_id",))
    record_json(name, f"{folder}/validation_summary.json",
                *both(folder, "validation_summary.json"))
    record_file(name, f"{folder}/BENCHMARK_REPRODUCTION.md",
                *both(folder, "BENCHMARK_REPRODUCTION.md"))
    record_file(name, "INTEGRATION_RUN.md", old / "INTEGRATION_RUN.md", new / "INTEGRATION_RUN.md")
    # the sequences, QC reports and per-base tables of the 13 samples, from the
    # evidence of the two regenerations
    before, after = args.prev_package / "evidence" / "mock", args.package / "evidence" / "mock"
    names = sorted({path.relative_to(root).as_posix()
                    for root in (before, after) for path in root.rglob("*")
                    if path.is_file() and path.name.endswith(MOCK_KINDS)})
    for relative in names:
        record_file(name, f"evidence/mock/{relative}", before / relative, after / relative)

    name, folder = "23S rRNA gene and gyrB", "low-copy-23S"
    table(name, f"{folder}/results.tsv", *both(folder, "results.tsv"), keys=("strain", "locus"))
    record_json(name, f"{folder}/evaluation.json", *both(folder, "evaluation.json"))
    record_json(name, f"{folder}/closed_truth_audit.json",
                *both(folder, "closed_truth_audit.json"))
    record_json(name, f"{folder}/run_summary.json", *both(folder, "run_summary.json"))
    record_json(name, f"{folder}/workflow_provenance.json",
                *both(folder, "workflow_provenance.json"))

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
    record_file(name, f"{folder}/README.md", *both(folder, "README.md"))

    name, folder = "LIBA-6656 consensus without a graph", "liba6656-depth"
    record_json(name, f"{folder}/consensus_without_graph.json",
                *both(folder, "consensus_without_graph.json"))
    table(name, f"{folder}/consensus_without_graph.tsv",
          *both(folder, "consensus_without_graph.tsv"))

    notes = []
    if args.notes and args.notes.is_file():
        for row in cr.read_table(args.notes):
            notes.append((re.compile(row["pattern"]), row["explanation"]))
    cases = rule_cases(args)
    rows = [row + [explain(row, notes, cases)] for row in comparison.rows]
    stem = f"old_vs_new_{args.previous_version}_vs_{args.version}"
    cr.write_table(args.package / f"{stem}.tsv", [
        "analysis", "record", "case", "quantity", "kind", f"value_{args.previous_version}",
        f"value_{args.version}", "equal", "absolute_difference", "explanation"], rows)
    counts = {}
    for row in rows:
        entry = counts.setdefault(row[0], {"values": 0, "equal": 0, "metadata": 0, "rule": 0,
                                           "explained": 0, "unexplained": 0})
        entry["values"] += 1
        if row[7] == "yes":
            entry["equal"] += 1
        elif row[9] == "UNEXPLAINED":
            entry["unexplained"] += 1
        elif row[9].startswith("run metadata"):
            entry["metadata"] += 1
        elif row[9].startswith("rule changed in 1.0.2"):
            entry["rule"] += 1
        else:
            entry["explained"] += 1
    cr.write_table(args.package / f"{stem}_summary.tsv", [
        "analysis", "values_compared", "equal", "different_run_metadata",
        "different_by_the_rule_changed_in_1.0.2", "different_explained_otherwise",
        "different_unexplained"],
        [[name, *entry.values()] for name, entry in counts.items()])
    for name, entry in counts.items():
        print(name, entry)
    print("RULE (all):")
    for row in [row for row in rows if row[9].startswith("rule changed")]:
        print("  ", " | ".join(str(value)[:60] for value in row[1:7]))
    print("NOT METADATA, NOT RULE (first 80):")
    for row in [row for row in rows if row[7] == "no" and not row[9].startswith(
            ("run metadata", "rule changed"))][:80]:
        print("  ", row[9][:14], "|", " | ".join(str(value)[:70] for value in row[:7]))


# -------------------------------------------------------------------- finish

def command_finish(args):
    """Keep in records/ only the files that differ from the 1.0.1 record."""
    records, old = args.package / "records", args.old
    diffs = args.package / "evidence" / "document_diffs"
    if diffs.exists():
        shutil.rmtree(diffs)
    rows = []
    for path in sorted(records.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(records)
        earlier = old / relative
        digest = cr.digest_file(path)
        if not earlier.is_file():
            rows.append([relative.as_posix(), "new file", "", digest, "yes"])
            continue
        before = cr.digest_file(earlier)
        if before == digest:
            rows.append([relative.as_posix(), "same bytes as the 1.0.1 record", before, digest, "no"])
            path.unlink()
            continue
        rows.append([relative.as_posix(), "differs from the 1.0.1 record", before, digest, "yes"])
        if path.suffix in (".md", ".txt"):
            diffs.mkdir(parents=True, exist_ok=True)
            lines = difflib.unified_diff(
                earlier.read_text().splitlines(), path.read_text().splitlines(),
                f"{args.previous_version}/{relative.as_posix()}",
                f"{args.version}/{relative.as_posix()}", lineterm="", n=0)
            (diffs / (relative.as_posix().replace("/", "__") + ".diff")).write_text(
                "\n".join(lines) + "\n")
    for directory in sorted((path for path in records.rglob("*") if path.is_dir()), reverse=True):
        if not any(directory.iterdir()):
            directory.rmdir()
    cr.write_table(args.package / "records_regenerated.tsv", [
        "file_under_validation", "status", f"sha256_{args.previous_version}",
        f"sha256_{args.version}", "in_records"], rows)
    for row in rows:
        print(row[4], row[1], row[0])


# ------------------------------------------------------------------- summary

def command_summary(args):
    package, old, records = args.package, args.old, args.package / "records"
    stem = f"old_vs_new_{args.previous_version}_vs_{args.version}"
    out = []

    def add(topic, quantity, value, source):
        out.append([topic, quantity, value, source])

    def record(relative):
        path = records / relative
        return path if path.is_file() else old / relative

    add("software", "Locus-Recon version", args.version, "PROVENANCE.json")
    add("software", "package source digest (SHA-256)", args.digest, "PROVENANCE.json")
    add("software", "git commit of the release tree", args.commit, "PROVENANCE.json")
    versions = [line.split("\t") for line in
                (args.regen / "COMMAND_VERSIONS.tsv").read_text().splitlines()]
    add("software", "version reported by the four commands",
        "; ".join(reported for _, reported in versions), "inputs/COMMAND_VERSIONS.tsv")
    files = cr.read_table(next(args.regen.glob("package_files_*.tsv")))
    differing = [row["file"] for row in files if row["same"] == "no"]
    add("software", f"files of locus_recon/ that differ between {args.previous_commit[:7]} and "
        f"{args.commit[:7]}", f"{len(differing)} of {len(files)} ({'; '.join(differing)})",
        f"inputs/{next(args.regen.glob('package_files_*.tsv')).name}")

    topic = f"comparison with the {args.previous_version} records"
    rows = cr.read_table(package / f"{stem}.tsv")
    old_column, new_column = f"value_{args.previous_version}", f"value_{args.version}"
    different = [row for row in rows if row["equal"] == "no"]
    by_rule = [row for row in different if row["explanation"].startswith("rule changed")]
    metadata = [row for row in different if row["explanation"].startswith("run metadata")]
    unexplained = [row for row in different if row["explanation"] == "UNEXPLAINED"]
    other = [row for row in different
             if row not in by_rule and row not in metadata and row not in unexplained]
    add(topic, "values compared", len(rows), f"{stem}.tsv")
    add(topic, "equal", len(rows) - len(different), f"{stem}.tsv")
    add(topic, "different run metadata", len(metadata), f"{stem}.tsv")
    add(topic, "different by the rule changed in 1.0.2", len(by_rule), f"{stem}.tsv")
    add(topic, "different explained otherwise", len(other), f"{stem}.tsv")
    add(topic, "different unexplained", len(unexplained), f"{stem}.tsv")
    for row in cr.read_table(package / f"{stem}_summary.tsv"):
        add(topic, f"{row['analysis']}: values compared / different",
            f"{row['values_compared']} / {int(row['values_compared']) - int(row['equal'])}",
            f"{stem}_summary.tsv")
    for row in by_rule:
        add(topic, f"changed by the rule: {row['case']}, {row['quantity']}",
            f"{row[old_column] or '(empty)'} -> {row[new_column] or '(empty)'}", f"{stem}.tsv")
    for row in other + unexplained:
        label = "unexplained" if row in unexplained else "explained otherwise"
        add(topic, f"{label}: {row['record']} | {row['case']} | {row['quantity']}"[:230],
            f"{(row[old_column] or '(empty)')[:70]} -> {(row[new_column] or '(empty)')[:70]}",
            f"{stem}.tsv")

    topic = "23S rRNA gene and gyrB"
    source = "records/low-copy-23S/results.tsv"
    earlier = {(row["strain"], row["locus"]): row
               for row in cr.read_table(old / "low-copy-23S" / "results.tsv")}
    for row in cr.read_table(record("low-copy-23S/results.tsv")):
        before = earlier[(row["strain"], row["locus"])]
        add(topic, f"{row['strain']} {row['locus']}: ratio (interval), depth call, graph count, "
            "consensus call, kind, method, status, lower bound, upper bound",
            f"{row['ratio']} ({row['ratio_ci_low']}-{row['ratio_ci_high']}), {row['depth_call']}, "
            f"{row['graph_context_count'] or '(none)'}, {row['copy_number_call'] or '(none)'}, "
            f"{row['copy_number_kind']}, {row['copy_number_method']}, {row['copy_number_status']}, "
            f"{row['copy_number_lower_bound'] or '(none)'}, "
            f"{row['copy_number_upper_bound'] or '(none)'}", source)
        add(topic, f"{row['strain']} {row['locus']}: status, lower bound, upper bound in "
            f"{args.previous_version}",
            f"{before['copy_number_status']}, {before['copy_number_lower_bound'] or '(none)'}, "
            f"{before['copy_number_upper_bound'] or '(none)'}",
            f"validation/low-copy-23S/results.tsv at {args.commit[:7]}")
    new_evaluation = json.loads(record("low-copy-23S/evaluation.json").read_text())
    old_evaluation = json.loads((old / "low-copy-23S" / "evaluation.json").read_text())
    for key, value in new_evaluation.items():
        add(topic, f"evaluation.json, {key}: {args.previous_version} -> {args.version}",
            f"{old_evaluation.get(key)} -> {value}", "records/low-copy-23S/evaluation.json")
    log = (args.regen / "logs" / "low_copy.log").read_text(errors="replace")
    found = re.findall(r"^RUN_LOW_COPY_23S_EXIT (\d+)", log, flags=re.M)
    add(topic, "exit code of the runner", found[-1] if found else "not found", "logs/low_copy.log")
    reused = cr.read_table(package / "evidence" / "copy_number" / "low_copy_23S" / "drafts_reused.tsv")
    add(topic, "draft files (contigs and graph) with the checksum used by the 1.0.1 run",
        f"{sum(row['identical'] == 'yes' for row in reused)} of {len(reused)}",
        "evidence/copy_number/low_copy_23S/drafts_reused.tsv")

    topic = "LIBA-6656 tcdB, standard workflow"
    source = "records/liba6656-gfa-rerun/RELEASE_VERIFICATION.json"
    verification = json.loads(record("liba6656-gfa-rerun/RELEASE_VERIFICATION.json").read_text())
    (report,) = cr.read_table(record("liba6656-gfa-rerun/standard_workflow_report_tcdB.tsv"))
    add(topic, "version in the verification record", verification["software"]["version"], source)
    add(topic, "length, tier, disposition, mixed sites",
        f"{report['allele_length']} bp, {report['sequence_confidence']}, "
        f"{report['result_disposition']}, {report['mixed_site_count']}",
        "records/liba6656-gfa-rerun/standard_workflow_report_tcdB.tsv")
    add(topic, "threads of the two runs; the runs agree",
        f"{' and '.join(str(run['parameters']['threads']) for run in verification['runs'])}; "
        f"{str(verification['runs_agree']).lower()}", source)
    add(topic, "graph equal to the deposited graph (segments, links)",
        f"{str(verification['graph']['same_segment_sequences']).lower()}, "
        f"{str(verification['graph']['same_links_with_segments_matched_by_sequence']).lower()}",
        source)
    add(topic, "terminal paths equal to the deposited paths",
        str(verification["terminal_paths"]["same_sequence_under_same_name"]).lower(), source)
    liba_rows = [row for row in rows if row["analysis"] == "LIBA-6656 standard workflow"]
    for label, pattern in (("reported sequence", "reconstructed_tcdB.fasta"),
                           ("QC report", "qc_report_tcdB.txt")):
        selected = [row for row in liba_rows if pattern in row["record"]]
        add(topic, f"{label} equal to the {args.previous_version} record",
            "yes" if selected and all(row["equal"] == "yes" for row in selected) else "no",
            f"{stem}.tsv")
    result = [row for row in liba_rows if not (
        row["equal"] == "yes" or row["explanation"].startswith("run metadata"))]
    add(topic, "values that differ and are not run metadata", len(result), f"{stem}.tsv")
    table = cr.read_table(package / "evidence" / "liba6656_1.0.2_vs_1.0.1_checksums.tsv")
    for label, suffix in (("local assembly graphs", "assembly_graph_after_simplification.gfa"),
                          ("reconstructed sequences", "_reconstructed.fasta"),
                          ("QC reports", "_qc_report.txt"),
                          ("per-base support tables", "allele_remap.support.tsv")):
        selected = [row for row in table if row["file"].endswith(suffix)]
        add(topic, f"{label} of the two runs with the bytes of the 1.0.1 runs",
            f"{sum(row['identical'] == 'yes' for row in selected)} of {len(selected)}",
            "evidence/liba6656_1.0.2_vs_1.0.1_checksums.tsv")

    topic = "LIBA-6656 tcdB, consensus without a graph"
    source = "records/liba6656-depth/consensus_without_graph.json"
    new_consensus = json.loads(record("liba6656-depth/consensus_without_graph.json").read_text())
    old_consensus = json.loads((old / "liba6656-depth" / "consensus_without_graph.json").read_text())
    for key in ("software_version", "copy_number_call", "copy_number_kind", "copy_number_method",
                "copy_number_status", "consensus_flags", "copy_number_lower_bound",
                "copy_number_upper_bound", "dosage_estimate", "dosage_status", "ratio",
                "depth_call", "geometry_expected_dosage"):
        add(topic, f"{key}: {args.previous_version} -> {args.version}",
            f"{old_consensus.get(key)} -> {new_consensus.get(key)}", source)
    selected = [row for row in rows if row["analysis"] == "LIBA-6656 consensus without a graph"]
    result = [row for row in selected if not (
        row["equal"] == "yes" or row["explanation"].startswith("run metadata"))]
    add(topic, "values compared / equal / run metadata / other",
        f"{len(selected)} / {sum(row['equal'] == 'yes' for row in selected)} / "
        f"{sum(row['explanation'].startswith('run metadata') for row in selected)} / {len(result)}",
        f"{stem}.tsv")

    topic = "simulated benchmark"
    summary = json.loads(record("reference_results/validation_summary.json").read_text())
    checked = cr.read_table(record("reference_results/validation_results.tsv"))
    report = cr.read_table(record("reference_results/locus_recon_report_mockLocus.tsv"))
    add(topic, "core criteria met",
        f"{sum(summary['criteria'].values())} of {len(summary['criteria'])}",
        "records/reference_results/validation_summary.json")
    add(topic, "samples", len(report), "records/reference_results/locus_recon_report_mockLocus.tsv")
    add(topic, "samples with a candidate identical to the truth",
        sum(row["exact_truth_match"].lower() == "true" for row in checked),
        "records/reference_results/validation_results.tsv")
    add(topic, "samples accepted (PASS)",
        sum(row["result_disposition"] == "PASS" for row in report),
        "records/reference_results/locus_recon_report_mockLocus.tsv")
    add(topic, "tool_version in the report", "; ".join(sorted({row["tool_version"] for row in report})),
        "records/reference_results/locus_recon_report_mockLocus.tsv")
    table = cr.read_table(package / "evidence" / "mock_1.0.2_vs_1.0.1_checksums.tsv")
    for label, suffix in (("reconstructed sequences", "_reconstructed.fasta"),
                          ("QC reports", "_qc_report.txt"),
                          ("per-base support tables", "allele_remap.support.tsv"),
                          ("pooled FASTA files", ("_alleles.fasta", "_candidates.fasta")),
                          ("local assembly graphs", "assembly_graph_after_simplification.gfa")):
        selected = [row for row in table if row["file"].endswith(suffix)]
        add(topic, f"{label} with the bytes of the {args.previous_version} run",
            f"{sum(row['identical'] == 'yes' for row in selected)} of {len(selected)}",
            "evidence/mock_1.0.2_vs_1.0.1_checksums.tsv")
    add(topic, f"all kept files with the bytes of the {args.previous_version} run",
        f"{sum(row['identical'] == 'yes' for row in table)} of {len(table)}",
        "evidence/mock_1.0.2_vs_1.0.1_checksums.tsv")
    repeat = cr.read_table(package / "evidence" / "mock_run_and_repeat_checksums.tsv")
    add(topic, "kept files with the same bytes in a second run of 1.0.2",
        f"{sum(row['identical'] == 'yes' for row in repeat)} of {len(repeat)}",
        "evidence/mock_run_and_repeat_checksums.tsv")
    selected = [row for row in repeat if row["file"].endswith(MOCK_KINDS[:3])]
    add(topic, "sequences, QC reports and support tables with the same bytes in the second run",
        f"{sum(row['identical'] == 'yes' for row in selected)} of {len(selected)}",
        "evidence/mock_run_and_repeat_checksums.tsv")

    topic = "unit tests"
    for path in sorted((package / "tests").glob("*.pytest.txt")):
        last = path.read_text(errors="replace").strip().splitlines()[-1]
        add(topic, path.name.replace(".pytest.txt", "").replace("_", " "), last.strip("= "),
            f"tests/{path.name}")
        failed = re.findall(r"^FAILED (\S+)", path.read_text(errors="replace"), flags=re.M)
        if failed:
            add(topic, path.name.replace(".pytest.txt", "").replace("_", " ") + ": failed tests",
                "; ".join(failed), f"tests/{path.name}")

    topic = "scope of the rerun"
    source = "inputs/reconcile_copy_number_reachability.tsv"
    reach = cr.read_table(args.regen / "reconcile_copy_number_reachability.tsv")
    column = "can_reach_reconcile_copy_number_by_imports_and_commands"
    left_out = [row for row in reach if row["rerun_with_1.0.2"] == "no"]
    add(topic, "analyses not rerun that reach reconcile_copy_number by imports and commands",
        f"{sum(row[column] == 'yes' for row in left_out)} of {len(left_out)}", source)
    add(topic, "analyses not rerun that reach it only if a script named in text is followed",
        "; ".join(row["analysis"] for row in left_out if row[column] == "no" and row[
            "further_referring_modules_if_scripts_named_in_text_are_followed"] != "none") or "none",
        source)
    add(topic, "analyses rerun that reach it",
        "; ".join(row["analysis"] for row in reach
                  if row["rerun_with_1.0.2"] == "yes" and row[column] == "yes"), source)
    add(topic, "analyses rerun that do not reach it",
        "; ".join(row["analysis"] for row in reach
                  if row["rerun_with_1.0.2"] == "yes" and row[column] == "no"), source)
    add(topic, "runners outside this task that reach it",
        "; ".join(row["analysis"] for row in reach
                  if row["rerun_with_1.0.2"] == "not in this task" and row[column] == "yes"),
        source)
    held = cr.read_table(args.regen / "records_not_rerun_that_hold_the_method.tsv")
    add(topic, "records not rerun that hold the method GRAPH_COUNT_OVER_DISCORDANT_DEPTH: "
        "files; lines; of these lines with SUPPORTED",
        f"{len(held)}; "
        f"{sum(int(row['lines_with_GRAPH_COUNT_OVER_DISCORDANT_DEPTH']) for row in held)}; "
        f"{sum(int(row['of_these_lines_with_SUPPORTED']) for row in held)}",
        "inputs/records_not_rerun_that_hold_the_method.tsv")
    add(topic, "folders of those records",
        "; ".join(sorted({"/".join(row["file"].split("/")[:2]) for row in held})) or "none",
        "inputs/records_not_rerun_that_hold_the_method.tsv")

    regenerated = cr.read_table(package / "records_regenerated.tsv")
    add("records", "record files regenerated / with the bytes of the 1.0.1 record / in records/",
        f"{len(regenerated)} / {sum(row['in_records'] == 'no' for row in regenerated)} / "
        f"{sum(row['in_records'] == 'yes' for row in regenerated)}", "records_regenerated.tsv")
    cr.write_table(package / f"summary_{args.version}.tsv",
                   ["topic", "quantity", "value", "source_file"], out)
    for row in out:
        print(" | ".join(str(value)[:150] for value in row[:3]))


# ------------------------------------------------------------------ manifest

def command_manifest(args):
    package, regen = args.package, args.regen
    tools = dict(line.split("\t", 1)
                 for line in (regen / "TOOL_VERSIONS.tsv").read_text().splitlines())
    commands_reported = dict(line.split("\t", 1)
                             for line in (regen / "COMMAND_VERSIONS.tsv").read_text().splitlines())
    steps, commands = [], {}
    for line in (regen / "logs" / "steps.tsv").read_text().splitlines():
        fields = line.split("\t")
        steps.append({"step": fields[0], **dict(field.split("=", 1) for field in fields[1:])})
    for log in sorted((regen / "logs").glob("*.log")):
        found = [line.lstrip("+ ").rstrip() for line in log.read_text(errors="replace").splitlines()
                 if re.match(r"^\++ (python|spades\.py|locus-recon|tar|git (clone|checkout|apply))", line)]
        if found:
            commands[log.stem] = found
    inputs = []
    for line in (regen / "inputs" / "INPUT_SHA256.tsv").read_text().splitlines():
        path, size, digest = line.split("\t")
        inputs.append({"path": path, "size_bytes": int(size), "sha256": digest,
                       "source": "ENA run" if path.endswith(".fastq.gz")
                       else "draft assembly of LIBA-6656 supplied by the author"})
    provenance = args.work / "low_copy_23S" / "results" / "workflow_provenance.json"
    if provenance.is_file():
        record = json.loads(provenance.read_text())
        for entry in record["inputs"]:
            inputs.append({**entry, "source": f"inputs recorded by {record['workflow']}"})
    status = subprocess.run(["git", "-C", str(args.repo), "status", "--porcelain"],
                            capture_output=True, text=True, check=True).stdout
    head = subprocess.run(["git", "-C", str(args.repo), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    cr.write_json(package / "PROVENANCE.json", {
        "schema_version": "1.0",
        "workflow": f"rerun with release {args.version} of the analyses that the release can change",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "software": {"name": "locus-recon", "version": args.version,
                     "source_sha256": args.digest, "git_commit": head,
                     "git_tree_clean": status == "",
                     "module_path": "the clone, placed first on PYTHONPATH; the conda "
                                    "environment holds an installed copy of release "
                                    f"{args.previous_version} that is not used"},
        "compared_with": {"version": args.previous_version, "git_commit": args.previous_commit,
                          "records": "validation/ of the release tree at the commit above"},
        "versions_reported_by_the_commands": commands_reported,
        "tools": tools, "steps": steps, "commands": commands,
        "seeds": {"bootstrap of the depth ratio (default of the commands and of the runners)":
                  20260817},
        "inputs": inputs,
        "reused_from_the_1.0.1_regeneration": [
            "work/tier/<strain>/spades/ (five SPAdes --isolate drafts, copied)",
            "work/low_copy_23S/closed_truth/ and fixed_baits/ (GenBank records, copied)",
            "work/liba6656_depth/results/ERR467623_vs_LIBA6656_draft.sorted.bam (symbolic link)",
            "inputs/LIBA6656_ST154.fasta (copied)", "pydeps/ (pytest, copied)",
        ],
        "driver": "scripts/validation/regenerate/regenerate_102.sh",
    })
    lines = []
    for path in sorted(package.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            lines.append(f"{cr.digest_file(path)}  {path.relative_to(package).as_posix()}")
    (package / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")
    total = sum(path.stat().st_size for path in package.rglob("*") if path.is_file())
    print(f"manifest: {len(lines)} files, {total} bytes; clone clean: {status == ''}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("evidence", "records", "integration", "compare",
                                            "finish", "summary", "manifest"))
    parser.add_argument("--regen", type=Path, required=True, help="working directory")
    parser.add_argument("--package", type=Path, required=True, help="package directory")
    parser.add_argument("--repo", type=Path, required=True, help="clone of the release")
    parser.add_argument("--previous", type=Path, required=True,
                        help="working directory of the 1.0.1 regeneration (read only)")
    parser.add_argument("--version", default="1.0.2")
    parser.add_argument("--previous-version", default="1.0.1")
    parser.add_argument("--digest", default="")
    parser.add_argument("--commit", default="")
    parser.add_argument("--previous-commit", default="")
    parser.add_argument("--notes", type=Path, default=None,
                        help="table of patterns and explanations of expected differences")
    args = parser.parse_args()
    args.work = args.regen / "work"
    args.old = args.repo / "validation"
    args.prev_package = args.previous / "package" / f"regenerated-{args.previous_version}"
    args.package.mkdir(parents=True, exist_ok=True)
    {"evidence": command_evidence, "records": command_records,
     "integration": command_integration, "compare": command_compare,
     "finish": command_finish, "summary": command_summary,
     "manifest": command_manifest}[args.command](args)


if __name__ == "__main__":
    main()
