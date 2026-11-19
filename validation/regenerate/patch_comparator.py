#!/usr/bin/env python3
"""Write the revised ``run_comparator_benchmark.py`` from the released one.

The released script took the rows of Locus-Recon from tables deposited under
``validation/``, because the candidate sequences of two data sets had not been
kept.  The revised script scores the candidate sequences of a regeneration
with the function used for ARIBA and KMA, and writes them beside the sequences
of the other tools.  Each change is one replacement of a passage that must
occur exactly once in the released file; the script stops otherwise.
"""

import argparse
from pathlib import Path

REPLACEMENTS = []


def replace(old, new):
    REPLACEMENTS.append((old, new))


replace(
    "``validation/``; Locus-Recon is not run here.\n",
    "``validation/``; Locus-Recon is not run here.\n"
    "\n"
    "With ``--locus-recon-evidence DIR`` the rows of Locus-Recon are not read from\n"
    "those tables.  They are scored from the candidate sequences kept by a\n"
    "regeneration of the validation analyses (``validation/regenerate``), against\n"
    "the same truth and with the same function as the sequences of ARIBA and KMA,\n"
    "and the candidates are written to ``comparator_returned_sequences.fasta``\n"
    "with those of the two tools.  ``DIR`` holds ``mock/`` (the results directory\n"
    "of the simulated benchmark), ``hpylori/<strain>/<locus>/`` (one output\n"
    "directory per closed-genome case) and ``liba6656/<run>/`` (the standard\n"
    "workflow on LIBA-6656).\n",
)

replace(
    "from validation.common import build_provenance, write_json  # noqa: E402\n",
    "from validation.common import build_provenance, sha256_file, write_json  # noqa: E402\n",
)

replace(
    "\n\ndef summarise(rows, groups):\n",
    '''

# ------------------------------------------- Locus-Recon rows scored from candidates

def _evidence_run(dataset, sample, locus, evidence):
    """Output directory of the Locus-Recon run of one case, and its sample name."""
    evidence = Path(evidence)
    if dataset == "hpylori":
        return evidence / "hpylori" / sample / locus, sample
    if dataset == "mock":
        return evidence / "mock", sample
    runs = sorted(path.parent for path in
                  (evidence / "liba6656").glob(f"*/locus_recon_report_{locus}.tsv"))
    if not runs:
        raise FileNotFoundError(f"no Locus-Recon run of {locus} under {evidence / 'liba6656'}")
    return runs[0], "LIBA6656_ST154"


def locus_recon_rows_from_candidates(dataset, samples, evidence, work):
    """Results of Locus-Recon for the same cases, scored from its candidates.

    Every candidate is compared with the truth by ``compare``, the function
    that scores the sequences of ARIBA and KMA, so that no field of a
    Locus-Recon row is derived from an earlier table.  Returns the rows, the
    candidate sequences and the files that were read.
    """
    rows, returned, files = [], [], []
    for sample in samples:
        for locus, truth in sample["truth"].items():
            run_dir, sample_id = _evidence_run(dataset, sample["sample"], locus, evidence)
            report = run_dir / f"locus_recon_report_{locus}.tsv"
            records = [record for record in calibration.read_tsv(report)
                       if record["sample_id"] == sample_id]
            if len(records) != 1:
                raise RuntimeError(f"{report}: expected one row for {sample_id}, "
                                   f"found {len(records)}")
            record = records[0]
            files.append(report)
            label = f"{record.get('tool_version') or 'unknown version'}, regenerated candidate"
            candidates = sorted((run_dir / sample_id).glob(f"*_{locus}_reconstructed.fasta"))
            name, sequence = "", ""
            if candidates:
                name, sequence = read_fasta(candidates[0])[0]
                files.append(candidates[0])
            row = {
                "dataset": dataset, "sample": sample["sample"], "locus": locus,
                "tool": LOCUS_RECON[0], "configuration": label,
                "reference_kind": REFERENCE_KIND[dataset], "truth_copies": len(truth),
                "sequences_returned": 1 if sequence else 0, "returned_bp": len(sequence),
                "uncalled_bases": sequence.upper().count("N"),
                "lowercase_bases": sum(1 for letter in sequence if letter.islower()),
                "truth_sequences_identical": recovered([sequence] if sequence else [], truth),
                "tool_reference": report.relative_to(Path(evidence)).as_posix(),
                "tool_report": _evidence_report(record, (
                    ("status", "workflow_status"), ("tier", "sequence_confidence"),
                    ("disposition", "result_disposition"), ("flags", "qc_flags"))),
            }
            if truth:
                row.update(truth_name=truth[0][0], truth_bp=len(truth[0][1]))
            if not sequence:
                row["comparison"] = "no sequence returned"
            elif truth:
                row.update(compare(sequence, truth, work))
                row["also_identical_to"] = ";".join(
                    other for other, other_seq in sample["other"].items()
                    if sequence.upper() in (other_seq, reverse_complement(other_seq)))
            else:
                row["comparison"] = RELATIONS[7]
            if sequence:
                identifier = (f"{dataset}|{sample['sample']}|{locus}|{LOCUS_RECON[0]}|"
                              f"{label.replace(', ', '_').replace(' ', '_')}|1")
                returned.append((identifier + " " + name, sequence))
            rows.append(row)
    return rows, returned, files


def summarise(rows, groups):
''',
)

replace(
    '''    parser.add_argument("--score-only", action="store_true",
                        help="use the raw outputs already under <workdir>/raw")
''',
    '''    parser.add_argument("--score-only", action="store_true",
                        help="use the raw outputs already under <workdir>/raw")
    parser.add_argument("--locus-recon-evidence", metavar="DIR",
                        help="score the Locus-Recon candidates kept under DIR "
                             "instead of reading the deposited tables")
''',
)

replace(
    "    rows, returned, truth_records, catalogues, removed, groups = [], [], [], [], {}, {}\n",
    "    rows, returned, truth_records, catalogues, removed, groups = [], [], [], [], {}, {}\n"
    "    evidence_files = []\n",
)

replace(
    '''        rows += locus_recon_rows(dataset, samples, workdir / "blast")
''',
    '''        if args.locus_recon_evidence:
            scored, sequences, files = locus_recon_rows_from_candidates(
                dataset, samples, args.locus_recon_evidence, workdir / "blast")
            rows += scored
            returned += sequences
            evidence_files += files
        else:
            rows += locus_recon_rows(dataset, samples, workdir / "blast")
''',
)

replace(
    "        tools = [LOCUS_RECON] + [(tool, name) for tool, name, _ in CONFIGURATIONS]\n",
    "        tools = sorted({(row[\"tool\"], row[\"configuration\"]) for row in rows\n"
    "                        if row[\"tool\"] == LOCUS_RECON[0]})\n"
    "        tools += [(tool, name) for tool, name, _ in CONFIGURATIONS]\n",
)

replace(
    '''    print(f"{len(rows)} rows written to {outdir}", flush=True)
''',
    '''    if args.locus_recon_evidence:
        evidence = Path(args.locus_recon_evidence).resolve()
        provenance = build_provenance(
            workflow="comparator_benchmark/locus_recon_rows_scored_from_candidates",
            repo_root=REPO_ROOT,
            inputs=[CALIBRATION_DIR / f"bait_{locus}_26695.fasta" for locus in HPYLORI_LOCI]
            + [CALIBRATION_DIR / "closed_genome_manifest.tsv",
               MOCK_DIR / "mock_locus_alleles.fasta", MOCK_DIR / "dataset_manifest.json",
               LIBA_BAITS, LIBA_CANDIDATES]
            + sorted(cache.glob("*_from_genomic.fna.gz")),
            parameters={"datasets": args.datasets.split(","), "score_only": args.score_only,
                        "locus_recon_evidence": evidence.name},
            tools={"python": [sys.executable, "--version"], "blastn": ["blastn", "-version"]},
        )
        provenance["locus_recon_files_read"] = [
            {"path": path.resolve().relative_to(evidence).as_posix(),
             "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in sorted(set(evidence_files))]
        write_json(outdir / "LOCUS_RECON_ROWS_PROVENANCE.json", provenance)
    print(f"{len(rows)} rows written to {outdir}", flush=True)
''',
)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--source", required=True, type=Path,
                        help="run_comparator_benchmark.py of the release")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    text = args.source.read_text()
    for old, new in REPLACEMENTS:
        count = text.count(old)
        if count != 1:
            raise SystemExit(f"passage found {count} time(s), expected once:\n{old}")
        text = text.replace(old, new)
    args.output.write_text(text)
    print(f"{len(REPLACEMENTS)} passages replaced; written to {args.output}")


if __name__ == "__main__":
    main()
