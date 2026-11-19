# Per-unit evidence

`SHA256SUMS` lists every file of the results folder as it was packed on the host. Two of its
subfolders are not kept in the repository because of their size:

- `per_unit/` (1,534 files): for every unit, the candidate FASTA, the QC report,
  `allele_remap.support.tsv`, the report row and the outputs of the two copy-number commands;
- `comparator_raw/`: the reports of ARIBA and KMA for every isolate.

Both are in Supplementary Data 1 of the article, folder
`validation_evidence/release_1.0.1/external-evaluation/`, with the same checksums.
