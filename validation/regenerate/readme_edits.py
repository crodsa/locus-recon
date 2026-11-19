"""Passages of README files that name the release or the source of a row.

``collect_regenerated.py records`` copies each README of the earlier release
and replaces these passages, each of which must occur exactly once.  The
numbers that the README files quote were all reproduced by the regeneration
(see ``old_vs_new.tsv``), so no number is changed here.
"""

EDITS = {
    "comparator-benchmark/README.md": [
        (
            "This folder holds a comparison of Locus-Recon 1.0.0 with two established tools "
            "that also return a locus\n"
            "sequence from short reads: ARIBA 2.14.7 and KMA 1.7.1. The three tools received "
            "the same reads and the same\n"
            "reference sequences, and the sequences they returned were scored against the "
            "same truth.\n",
            "This folder holds a comparison of Locus-Recon with two established tools "
            "that also return a locus\n"
            "sequence from short reads: ARIBA 2.14.7 and KMA 1.7.1. The three tools received "
            "the same reads and the same\n"
            "reference sequences, and the sequences they returned were scored against the "
            "same truth. The sequences of\n"
            "ARIBA and KMA are those of their runs of 11 and 13 November 2026. The sequences "
            "of Locus-Recon are those of\n"
            "the regeneration of the validation analyses with release 1.0.1, of 18 November "
            "2026.\n",
        ),
        (
            "It runs ARIBA and KMA. The rows\n"
            "of Locus-Recon are read from the tables that are already in this repository; "
            "Locus-Recon is not run again.\n",
            "It runs ARIBA and KMA. With the\n"
            "option `--locus-recon-evidence` it scores the candidate sequences that the "
            "regeneration kept\n"
            "(`../regenerated-1.0.1/evidence/`) with the function that scores the sequences "
            "of the other two tools.\n"
            "Locus-Recon itself is not run by this script. Without that option the rows of "
            "Locus-Recon are read from\n"
            "tables of this repository, as in release 1.0.0.\n",
        ),
        (
            "  13 November 2026 (ARIBA on the simulated samples), from source distributions "
            "of release 1.0.0 with the same\n"
            "  package code. The scoring and `PROVENANCE.json` date from 13 November 2026.\n",
            "  13 November 2026 (ARIBA on the simulated samples), from source distributions "
            "of release 1.0.0 with the same\n"
            "  package code. `PROVENANCE.json` and the first scoring date from 13 November "
            "2026.\n"
            "- On 18 November 2026 the outputs of ARIBA and KMA kept in "
            "`raw_tool_outputs.tar.gz` were scored again, without\n"
            "  running the two tools, together with the Locus-Recon candidates of release "
            "1.0.1\n"
            "  (`LOCUS_RECON_ROWS_PROVENANCE.json`). Every field of the ARIBA and KMA rows "
            "was reproduced.\n",
        ),
        (
            "The Locus-Recon\n"
            "candidates of the *H. pylori* and simulated cases are not stored in this "
            "repository, so their relation to the\n"
            "truth is derived from the deposited comparison fields.\n",
            "In release 1.0.0 the Locus-Recon\n"
            "candidates of the *H. pylori* and simulated cases had not been kept, and their "
            "relation to the truth was\n"
            "derived from deposited comparison fields. They are now scored from the "
            "sequences, which are in\n"
            "`comparator_returned_sequences.fasta`. `../regenerated-1.0.1/old_vs_new.tsv` "
            "lists every field of the earlier\n"
            "and of the new rows.\n",
        ),
        (
            "| `comparator_returned_sequences.fasta` | Every sequence returned by ARIBA and "
            "KMA |\n",
            "| `comparator_returned_sequences.fasta` | Every sequence returned by ARIBA and "
            "KMA, and every candidate of Locus-Recon |\n"
            "| `LOCUS_RECON_ROWS_PROVENANCE.json` | The Locus-Recon reports and candidates "
            "that were scored, with their checksums |\n",
        ),
        (
            "python validation/run_comparator_benchmark.py --workdir work --outdir work/out "
            "--score-only\n```",
            "python validation/run_comparator_benchmark.py --workdir work --outdir work/out "
            "--score-only\n"
            "\n"
            "# the same, with the Locus-Recon rows scored from the candidates of the "
            "regeneration\n"
            "python validation/run_comparator_benchmark.py --workdir work --outdir work/out "
            "--score-only \\\n"
            "    --locus-recon-evidence validation/regenerated-1.0.1/evidence\n```",
        ),
    ],
}
