# Depth-based copy number: prespecified criteria and validation record

> **Editorial note.** Sections 1 to 5 record what was fixed in advance of the
> analyses; section 6 is post hoc and says so. Where Locus-Recon 1.0.0 departs
> from the prespecification, the departure is listed here rather than written
> into the record:
>
> - The windowed-profile criteria (dispersion index H, window P90 and the
>   4 kb/20-window calibration geometry) are not implemented. They were
>   calibrated against single-copy backbone segments and never evaluated
>   against loci with known partial amplifications, so their sensitivity was
>   never measured.
> - Incomplete coverage of the query does not withhold `dosage_estimate`: the
>   estimate is reported as the mean multiplicity over the covered query bases,
>   with `dosage_status` `ESTIMATED_PARTIAL_COVERAGE`, the covered fraction and
>   a bracket in which the uncovered bases count as zero copies or as the
>   largest multiplicity observed.
> - The external-validation record in section 4 (`external_validation.json`) is
>   the run that fixed this record. The released estimator resamples at most
>   512 evenly spaced backbone blocks for its interval; its results for the same
>   seven runs are in
>   [`../aphA1-copy-number/results.tsv`](../aphA1-copy-number/results.tsv).
>   Point estimates and calls are identical (`single-copy consistent` is
>   reported as `SINGLE_COPY_COMPATIBLE`), and the released intervals are
>   slightly wider, by at most 1.2 copies at 78 copies.
>
> The aggregate ratio criteria and the external qPCR comparison are in force as
> written and are what the released code applies.


This directory documents how the decision rules in `locus_recon/depth_ratio.py`
were fixed, what null distributions their thresholds come from, and how the
module behaves on external data with independently published copy numbers.

Every threshold below was set from a null distribution or from a written
criterion *before* the corresponding data were analysed. No threshold was
adjusted after seeing a result.

## 1. What the module decides, and how

The module compares filtered depth over the discovery intervals of a locus with
filtered depth over a single-copy backbone:

    ratio = median depth(locus) / median depth(backbone)

Numerator and denominator use identical filters (base quality >= 20, mapping
quality >= 20, `UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY` excluded). Contigs
shorter than 5 kb and 100 bp at every contig end are kept out of the backbone,
as are any regions passed to `exclude_regions` (rRNA operons, IS elements,
plasmid replicons — a replicon with its own copy number would bias the
denominator). The confidence interval is a block bootstrap over contiguous
500 bp blocks, 1,000 resamples and a fixed seed, because per-base depth is
autocorrelated and resampling single positions would give an absurdly tight
interval. A naive implementation of that bootstrap sorts the complete
multi-megabase backbone once per resample, which is not executable at
chromosome scale. The point denominator uses every eligible base. CI resampling
is capped at 512 backbone blocks
selected at deterministic, evenly spaced indices; the total and used block
counts are written to every result. This computational cap is fixed across
datasets and does not alter either decision threshold.

The aggregate and local-profile decisions are independent:

| estimand | condition | call |
|---|---|---|
| aggregate depth | `ratio_ci_low > call_threshold` (default 1.5) | `MULTICOPY_DEPTH` |
| aggregate depth | otherwise | `SINGLE_COPY_COMPATIBLE` |
| calibrated profile | window `P90 > p90_threshold` (default 1.65) | `LOCAL_EXCESS` |
| calibrated profile | otherwise | `NO_LOCAL_EXCESS` |
| any other profile geometry | descriptive values only | `NOT_EVALUATED_UNCALIBRATED_GEOMETRY` |

`ambiguity_index` (filtered / unfiltered depth at the locus) annotates the
result but never blocks it. The bias it measures runs one way: filtered depth
can never exceed raw depth, and a high-MAPQ single-copy backbone is not
depressed by the filter, so a low ambiguity index can only make the ratio too
small. Using it to withhold a rejection of single copy would discard precisely
the conclusion the bias protects. Below `ambiguity_warn` (0.70) the ratio is
therefore reported as a **lower bound**.

The dispersion index `H = (P75 - P25) / P50` over window ratios is reported as
a diagnostic and carries no verdict. It responds to dispersion of any origin —
MAPQ collapse over repeated flanks and the genuine edges of an amplification
alike — and is blind to sign: in the *aphA1* series below it reaches 0.644 on a
single-copy control and 2.104 on a high-copy amplification, so no threshold on H
separates the two classes.

## 2. Null distributions behind the window thresholds

Single-copy backbone segments of 4,000 bp (20 windows of 200 bp), quality
filters as above, taken from a genome unrelated to the validation panel:

| statistic | n segments | P50 | P90 | P95 | P99 | max |
|---|---|---|---|---|---|---|
| H | 916 | — | — | — | **0.288** | — |
| window P90 | 916 | — | — | — | **1.617** | — |

Thresholds were rounded up from the P99 of the null to 0.30 and 1.65. Both were
derived from single-copy sequence only; no amplified locus contributed to them.
`scripts/null_distribution.py` regenerates both from a BAM.

Because the calibration scale is 20 windows, a shorter single-region locus is
profiled over a 4 kb interval centred on the gene only when the BAM header
proves that the span remains inside the same contig. Users may instead provide
an explicit, independent `profile_region`. Calibration requires exactly one
4,000 bp interval and twenty complete, eligible 200 bp windows. A missing,
dropped, partial, shorter, longer, or multi-region profile remains reportable
but cannot use the fixed P90 decision. The aggregate ratio, its interval and
the ambiguity index are always computed on `locus_regions`, independently of
the profile.

## 3. Length-weighted copy estimate

The ratio measures **collapse**: the pile-up of reads onto a representative the
assembler failed to separate. Where the assembler resolved the copies onto
different contigs there is no pile-up, and a ratio near 1 is the correct answer
to the question the ratio asks while saying nothing about how many copies exist.
`dosage_estimate` combines resolved length with collapsed depth:

    copies = SUM_r (aligned_gene_bp_r * ratio_r) / gene_length

where `aligned_gene_bp_r` comes from the **query** coordinates of the discovery
BLAST — how much of the bait each assembly fragment covers — not from the padded
subject intervals.

Four predictions were written down before any of them was evaluated.
`scripts/estimator_predictions.py` reproduces all four through the released
function itself, so the number quoted below and the number a reader obtains
come from the same code path; the machine-readable output is
`estimator_predictions.json`.

| # | prediction | required | observed | verdict |
|---|---|---|---|---|
| a | on a single-fragment locus the estimator reduces exactly to `ratio` | deviation 0.00 on all 5 single-fragment *aphA1* cases | 0.95, 1.09, 1.10, 10.45, 77.67; deviation 0.00 | PASS |
| b | on single-copy 4 kb backbone segments the estimator stays at 1 | median in [0.90, 1.10], P99 <= 1.65 | n = 925; median 0.984, P99 1.226, max 1.274 | PASS |
| e | splitting a single-copy locus into n fragments must not multiply the estimate | slope over n below 0.02 copies/fragment | n = 512 stretches of 7,104 bp split into n = 2..9: 0.984–0.989, slope **-0.00022** | PASS |
| f | two disjoint single-copy stretches must give 2 | median near 2 | n = 256; median 1.985, P25/P75 1.81/2.17 | PASS |

Test (e) is the one that matters. The naive alternative
`ratio * assembly_copies` uses a misleadingly named fragment count and returns
1.97 (n = 2) through 8.85 (n = 9) on the same data, where
the truth is 1 copy: it counts assembly fragments, not copies. When a locus
spans several regions and no gene coordinates are supplied, `dosage_estimate`
is left at `None` rather than guessed. The same applies when the union of query
spans does not cover the full query. Overlapping query spans are not silently
deduplicated: total aligned, union-covered, and overlapping bases are reported;
the additive estimate is retained with `REVIEW_OVERLAPPING_SPANS`.
`assembly_copies` and `copies_estimate` carry the same values as
`assembly_region_count` and `dosage_estimate`.

## 4. External validation: *aphA1* in *Acinetobacter baumannii*

Seven public Illumina runs from strains whose *aphA1* copy number was determined
independently by qPCR and by WGS coverage. Reference `GCA_019458485.1`
(MRSN 56); bait `CP080452.1:250422-251237`, 816 bp, single BLAST hit at 100%
identity; a co-annotated *APH(6)-I* hit was rejected as a different gene. The
four plasmid replicons were excluded from the backbone. Expected classes were
fixed in a manifest before any run was analysed. Results in
`external_validation.json`; read accessions and MD5 sums in
`download_provenance.tsv`.

| strain | run | expected | ratio (95% CI) | call | independent truth |
|---|---|---|---|---|---|
| MRSN 56 | SRR14998418 | ~1x | 1.098 (1.088–1.104) | single-copy consistent | 1 copy, MIC 0.5 |
| MRSN 3364 | SRR15734211 | ~1x | 1.017 (0.959–1.048) | single-copy consistent | 1 copy, MIC 0.5 |
| MRSN 3363 | SRR15734212 | ~1x | 1.087 (1.038–1.115) | single-copy consistent | 1 copy, MIC 0.5 |
| MRSN 3361 | SRR15734215 | ~1x | 0.952 (0.945–0.952) | single-copy consistent | 1 copy, MIC 0.5 |
| MRSN 57 | SRR17023755 | >1x | 10.451 (10.211–10.771) | MULTICOPY_DEPTH | 10 +/- 2 by qPCR, MIC 8 |
| MRSN 58 | SRR17023784 | >1x | 77.672 (74.180–78.459) | MULTICOPY_DEPTH | 75 +/- 14 by qPCR, MIC 16 |
| MRSN 56 day 10 | SRR15115349 | >1x | 91.408 (88.799–96.154) | MULTICOPY_DEPTH | amplified under selection; no published point value |

Class concordance was 4/4 among single-copy runs and 3/3 among amplified runs.
The aggregate gene-depth estimates for MRSN 57 and MRSN 58 fall within their
corresponding published qPCR uncertainty intervals. The gene interval is the
aggregate estimand; the separate 4 kb interval was the only geometry eligible
for the prespecified profile call, which the released software does not make.

## 5. Declared limitations

- The module separates one copy from two with margin; it does not claim to
  separate two from 2.5. The resolution between adjacent low copy numbers is not
  established by this evidence; how far it extends was measured afterwards and
  is reported in section 6.
- No external, independently published two-copy truth set was available, so the
  two-copy regime is supported by the constructed geometries of test (f) and by
  a single biological case with a closed-genome reference.
- The ratio is a mean over the sequenced population. An amplification carried by
  a subpopulation gives a fractional value, which is reported as such and not
  rounded.
- Copy number of a plasmid-borne gene is not inferred from the depth of that
  gene alone.

## 6. Post hoc resolution between adjacent copy numbers

This section is post hoc. It was written after the results of section 3 were
known, answers a different question, and could not have falsified anything.
No criterion, threshold or observed value above was changed by it.

The question is how far a point estimate near two copies constrains the
underlying count. Geometries carrying exactly one, two and two and a half
copies of a 7,104 bp gene were cut out of single-copy backbone sequence, the
released `length_weighted_copies` was applied to each, and the estimate of
2.258 obtained on the *tcdB* locus was placed inside the three resulting
distributions as a two-sided rank. Full output in
`copy_number_resolution.json`.

| truth | n geometries | median | P25 | P75 | two-sided p for 2.258 |
|---|---|---|---|---|---|
| 1.0 | 512 | 0.988 | 0.899 | 1.081 | 0.0039 |
| 2.0 | 256 | 1.985 | 1.807 | 2.169 | 0.180 |
| 2.5 | 170 | 2.466 | 2.257 | 2.721 | 0.506 |

One copy is excluded. Two is not, and neither is two and a half: the estimate
sits almost exactly at the median of the 2.5-copy distribution, so it is at
least as compatible with 2.5 as with 2. A midpoint rule at 2.25 would
misassign 10.5% of true two-copy geometries and 22.4% of true 2.5-copy ones,
which is why no such rule is applied and why the limitation above stands as
written.

The same construction bounds the upper tail at one copy. The largest estimate
obtained from single-copy sequence was 4.02, on contig `.11940_5_90.3` at
position 1, where one of the four fragments falls on a backbone region at
12.85 times the median depth. Three of the four fragments sit at 1.03 to 1.10,
so the aggregate ratio stays near 1 and the geometry is not called
`MULTICOPY_DEPTH`. H = 2.075 and window P90 = 12.94 expose the local
heterogeneity, but a fragmented multi-region profile is outside the fixed
4 kb/20-window calibration and therefore receives no P90 verdict. This case
motivates reporting geometry and profile calibration state explicitly rather
than allowing a profile statistic to overwrite the aggregate depth call.

## 7. Files

| file | content |
|---|---|
| `external_validation.json` | full module output for the seven *A. baumannii* runs, both intervals |
| `download_provenance.tsv` | ENA accessions, file names and MD5 sums of the reads used |
| `estimator_predictions.json` | machine-readable results of predictions (a), (b), (e), (f) |
| `copy_number_resolution.json` | post hoc resolution distributions at one, two and 2.5 copies |
| `estimator_end_to_end.json` | length-weighted estimate on a locus split across nine assembly fragments |
| `scripts/null_distribution.py` | regenerates the H and P90 null distributions from a BAM |
| `scripts/estimator_predictions.py` | evaluates predictions (a), (b), (e), (f) through the released function |
| `scripts/copy_number_resolution.py` | rebuilds the post hoc resolution distributions of section 6 |

Field names in the deposited JSON files were translated into English; the values
are the unmodified output of the analysis runs.
