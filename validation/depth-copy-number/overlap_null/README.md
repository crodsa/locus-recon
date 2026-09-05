# Overlap null control for the additive dosage

## Question

The additive per-query-base dosage sums each discovery interval's ratio over the
query bases it resolves. When intervals overlap, a query base is resolved several
times and accumulates the sum of their ratios. Does the resulting number measure
copy number, or does it measure the overlap?

## Design

`run_overlap_null.py`. Pseudo-queries of 7,104 bp — the length of the LIBA-6656
*tcdB* query — are drawn from repeat-rich positions of the same draft assembly,
located by self-alignment. The eight contigs carrying the locus are excluded.
Each pseudo-query is passed through the same discovery alignment and the same
dosage estimand as a real query, using the same BAM, filters and seed. Windows
returning a single discovery interval are discarded, since they cannot test
overlap. Of 3,000 windows tried, 87 returned two or more intervals.

None of the retained queries has independent evidence of duplication. They are
stretches of a draft assembly. Any dosage above one that they return is produced
by the geometry and the mapping, not by copy number.

## Result

Dosage rises with the fraction of the query resolved by more than one interval:

    dosage = 0.930 + 1.539 * overlap_fraction     R^2 = 0.636, residual SD 0.135

Observed overlap fractions span 0.049 to 0.453. The *tcdB* query's 0.833 is
outside that range, so the null at that geometry is an extrapolation: 2.211 with
a 95% prediction interval of [1.894, 2.528], against an observed 2.287 (z = 0.48).

The purely combinatorial expectation — the same spans with every region at unit
ratio — is 1.843. The empirical control lands higher because multiply represented
sequence also attracts multi-mapped reads, which raises the per-region ratios in
addition to the combinatorial effect. The combinatorial figure is therefore a
lower bound on the inflation, and it is what `geometry_expected_dosage` reports.

## Consequence for the software

`dosage_status` returns `REVIEW_OVERLAPPING_SPANS` whenever spans overlap. No
correction is applied and no threshold is fitted: this control comes from one
genome, one depth and one aligner, and calibrating a correction on a single
geometry is the failure mode that removed the windowed depth profile from this
release. The control establishes that the estimand is not identified under
overlap; it does not establish a transferable slope.

## Files

- `run_overlap_null.py` — the control
- `null_dose_overlap_control.csv` — 87 retained pseudo-queries
- `null_dose_overlap_control.json` — fit, prediction and the case comparison
