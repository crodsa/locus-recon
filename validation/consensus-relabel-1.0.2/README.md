# Consensus labels of release 1.0.2 applied to the results of release 1.0.1

The external evaluation (`validation/external-evaluation/`) and the targeted tests
(`validation/targeted-tests/`) were run with release 1.0.1, as their protocol states, and their
records are kept as they were written. The external evaluation showed that release 1.0.1 gave the
status `SUPPORTED` to integer graph counts that the depth interval excludes, and that five of
these counts were wrong (rRNA genes with six or seven copies and two contexts in the graph).
Release 1.0.2 returns such a count with the status `REVIEW`, as a lower bound and without an
upper bound. No other rule changed.

`validation/relabel_consensus.py` reads the stored depth and graph result of every unit, calls
`reconcile_copy_number` of the installed release on them and writes the stored consensus beside
the recomputed one.

| File | Content |
|---|---|
| `consensus_relabel.tsv` | One row per unit and configuration (678 rows): 128 units of the external evaluation, and the 275 samples of targeted tests 4 and 5 with and without the graph. Stored consensus (1.0.1), recomputed consensus (1.0.2), the fields that differ, each scored against the true copy number, and whether the recomputed lower bound holds |
| `consensus_relabel_summary.tsv` | Units by analysis, status and score, for the stored and the recomputed consensus |
| `PROVENANCE.json` | Version and digest of the package that recomputed the consensus; checksums of the three input tables |
| `check_with_1.0.1.json` | The same script run with release 1.0.1 (tag `v1.0.1`) and `--expect-identical`: 678 of 678 stored results reproduced. This shows that the stored fields are sufficient to recompute the consensus |

## Result

Sixteen of 678 results change, all in `copy_number_status` (`SUPPORTED` to `REVIEW`) and
`copy_number_upper_bound` (the count to empty). No call, kind, method or flag changes.

- **External evaluation, protocol set (113 units with a truth copy).** With 1.0.1, 83 of 89
  `SUPPORTED` results were correct, 6 of 14 `REVIEW` and 10 of 10 `LOWER_BOUND`. With 1.0.2,
  81 of 82 `SUPPORTED` results are correct, 8 of 21 `REVIEW` and 10 of 10 `LOWER_BOUND`. Seven
  results move to `REVIEW`: the five wrong counts of two (true copy numbers six and seven) and
  two correct counts of three (*Bordetella pertussis* FDAARGOS_1365, 23S and 16S rRNA genes,
  depth ratios 4.47 and 4.38). All seven counts hold as lower bounds. The `SUPPORTED` result that
  remains wrong is the 16S rRNA gene of FDAARGOS_209, the isolate whose reads do not match its
  genome.
- **Targeted tests 4 and 5 (275 samples).** Nine correct counts move to `REVIEW` (seven in
  test 4, two in test 5), because the bootstrap interval is narrower than the error of the
  ratio. No result was `SUPPORTED` and wrong before or after.

The rule therefore removes the wrong `SUPPORTED` counts of the evaluation at the cost of 11
correct counts that are now returned for review. It was written after these results were known
and has not been tested on other data.

The four 23S rRNA counts of *Helicobacter pylori* in `validation/low-copy-23S/` were run again
with release 1.0.2 and are `REVIEW` in that record.
