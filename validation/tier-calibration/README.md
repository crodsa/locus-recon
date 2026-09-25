# Tier calibration against known truth

Every other validation stage asks whether a reconstruction is correct. This one
asks what a reported tier buys the reader, because that is the question a
reviewer asks of any tool that withholds results: if the top tier is rarely
reached, is it reachable at all, and when it is reached, is it right?

Two quantities answer it, and only the first is a property of the tool.

- **Precision of the top tier**: of the results reported as deposit-ready, how
  many match truth exactly. One incorrect result in that tier would make the
  tier worthless, so the count of false accepts is the primary number.
- **Acceptance rate**: what fraction of results reach the top tier. This is a
  property of the case mix. A benchmark built mostly from degraded inputs is
  expected to accept few of them, so the rate is reported per case class, and
  the classes are prespecified from the input, never read off the results.

## Case classes

| class | definition | correct behaviour |
|---|---|---|
| `supported` | single-copy locus within 95% identity of the bait allele, intact in the draft, adequate depth, pure culture | accept |
| `divergent` | single-copy locus further than 95% identity from the one-allele bait | accept: the tier reports read support for the reported bases, not similarity to a reference |
| `low-depth` | the same single-copy locus reconstructed from reads subsampled to a fraction of their depth | tier falls with depth; a degraded reconstruction is withheld, not accepted |
| `truncated` | the locus runs off a contig end by construction | tier falls with the size of the loss |
| `multi-copy` | the closed genome annotates more than one copy, so near-identical copies collapse in a short-read draft | withhold, even when the reported consensus is exact |
| `degraded` | depth below the stated floor, a controlled allele mixture, a paralogue, or a target-negative genome | withhold |

A locus is assigned to `divergent` by its measured identity to the bait, not by
which arm produced it: the selection step ranks candidates by distance, and two
of the three loci it returned turned out to sit as close to the bait as the
housekeeping locus does, so they are counted as `supported`.

## Result

62 cases with known truth, from four stages, 35 of them on real reads. The
per-class tables of every stage are rendered from `TIER_CALIBRATION.json` into
`tier_calibration_tables.md` by the run itself.

### All four stages

[[CALIB-13 per-class table for all four stages, pasted from tier_calibration_tables.md]]

Four statements answer the question the deposit asks. The first three are
settled on the deterministic stage below and are tested again on real reads
here; the fourth needs the real-read stages.

1. **Is the top tier reachable, and is it precise?** Supported cases that
   reached `HIGH`: [[CALIB-05]]. Precision of `HIGH`: [[CALIB-06]]. False
   accepts: [[CALIB-07]].
2. **An overall acceptance rate is not a performance measure.** Most cases
   were built or subsampled to be refused, so the rate is a property of the
   case mix. On a dataset of fragmented or hypervariable targets it is expected
   to approach zero, and that is the designed behaviour rather than a failure.
3. **A withheld result is not a wrong result.** Withheld cases that
   reconstructed truth exactly over the span they reported: [[CALIB-08]]. The
   tier states what the reads establish, not what the sequence happens to be,
   so a withheld result means "not established", never "incorrect".
4. **Where the reported bases disagree with truth, what tier did they get?**
   [[CALIB-09]]. For every disagreement, `relation_to_annotation` in the
   per-case tables separates a base difference from a boundary difference.

## Stage A: deposited benchmarks (27 cases, deterministic)

Recomputed from `../reference_results/validation_results.tsv` (13-sample mock
truth set) and `../completeness-benchmark/completeness_benchmark_results.tsv`
(14 constructed completeness cases), both of which carry a per-case truth
verdict alongside the tier. Per-case table:
`tier_calibration_benchmarks.tsv`.

| case class | reached the top tier | top-tier precision (95% CI) | false accepts | withheld but exact |
|---|---|---|---|---|
| `supported` | 6 / 6 | 1.00 (0.54-1.00) | 0 | 0 |
| `truncated` | 0 / 11 | — | 0 | 11 |
| `degraded` | 0 / 10 | — | 0 | 8 |
| **all** | 6 / 27 | 1.00 (0.54-1.00) | 0 | 19 |

Six of six `supported` cases reached `HIGH` with exact truth, including the two
where the locus was broken across contigs and rejoined. All eleven `truncated`
cases were withheld and reconstructed exactly over the recovered span: the
eight simple truncations and the split whose missing interval the assembly
lacks at `MEDIUM`, the two splits whose missing interval another contig
carries at `SUSPECT`. The class table expected the tier to fall with the size
of the loss; it does not fall further once a loss is reported. Any reported
loss blocks `HIGH` and is scored once, and only a loss whose sequence another
contig carries is held. Eight of the ten `degraded` cases also returned the
exact allele while being withheld, among them the 3x and 5x mock samples at
`SUSPECT`. The paralogue control and the target-negative genome are the only
cases here that do not match truth, and both were withheld.

## Stage B: closed-genome audit on real reads (10 cases)

Five closed *Helicobacter pylori* chromosomes with their public Illumina runs.
The bait is a fixed external catalogue from strain 26695; sample-specific truth
is never used as a bait. Each draft is assembled with SPAdes `--isolate`, each
locus reconstructed, and each candidate compared to the annotated gene of its
own closed genome. The 23S rRNA gene is reconstructed with `--noncoding-locus`,
so reading-frame checks are not applied to it. Accessions:
`closed_genome_manifest.tsv`; per-case results:
`tier_calibration_closed_genomes.tsv`.

| Quantity | Result |
|---|---|
| *gyrB* loci at `HIGH` that match the annotated gene exactly | [[CALIB-02]] of 5 |
| Identity of the *gyrB* reconstructions to the nearest catalogue allele | [[CALIB-03]] |
| Two-copy 23S loci: tiers, and consensus exact to the annotation | [[CALIB-12]] |
| Reconstructions that were not exact, and how they were reported | [[CALIB-04]] |

## Stage C: divergent single-copy loci (15 cases)

The same drafts, with loci chosen by measurement rather than by hand. The
shortlist is prespecified in `CANDIDATE_LOCI`: single-copy *H. pylori* genes
with no known paralogue family, spanning housekeeping conservation to the
mosaic virulence genes. A locus is eligible only if the bait assembly and all
five closed genomes annotate it exactly once, and the three eligible loci
furthest from the bait allele are used. Every candidate and its per-genome copy
count is deposited in `locus_selection.tsv`; per-case results in
`tier_calibration_divergent_loci.tsv`.

[[CALIB-11 per-locus table: mean identity to the bait, bait length, cases at HIGH, relation to truth]]

A locus is counted as `divergent` only when its measured identity to the bait
is below 95%; a selected locus that sits as close to the bait as the
housekeeping locus is counted as `supported`. Where a divergent locus is
withheld, `relation_to_annotation` separates a response to wrong bases from a
response to an uncertain boundary.

Loci in the shortlist that the bait assembly or the closed genomes do not
annotate under the same gene name are dropped for lack of a comparable
annotation, not for biology; the counts are in `locus_selection.tsv`.

## Stage D: depth series (10 cases)

*gyrB*, the locus accepted at full depth, reconstructed again from reads
subsampled to ~15x and ~8x. Subsampling keeps every k-th pair, which is
deterministic and preserves the pairing, and the draft is reassembled at each
depth so the loss is felt by the assembly as well as by the reconstruction.
Per-case results: `tier_calibration_depth_series.tsv`.

[[CALIB-10 per-depth table: tiers and exact sequences at full depth, ~15x and ~8x]]

The question is whether the response is graded and monotone: whether a library
that lost depth kept its tier or was accepted, and whether a reported sequence
with a wrong base, if any, was withheld.

**The arms differ in assembler mode.** The full-depth drafts use SPAdes
`--isolate`; the subsampled drafts use SPAdes' default mode, because `--isolate`
is documented for high-coverage isolate data. Each arm therefore uses the
settings one would actually use at that depth, and the comparison across
depths includes that difference. The control below separates the two.

## Command

```bash
# Stage A only: deterministic, no network or external tools
python validation/run_tier_calibration.py --outdir validation/tier-calibration

# All four stages and the assembler-mode control
python validation/run_tier_calibration.py \
  --outdir validation/tier-calibration \
  --manifest validation/tier-calibration/closed_genome_manifest.tsv \
  --workdir /scratch/tier-calibration --threads 24 \
  --divergent-loci 3 --depth-series 15,8 --depth-isolate
```

## Assembler-mode control

The same subsampled libraries are reassembled with `--isolate` at every depth
and `gyrB` reconstructed again (`tier_calibration_depth_series_isolate.tsv`).
Tier agreement with the primary series:
[[CALIB-14 tier agreement of the --isolate control with the primary depth series]].
The control reuses the same libraries and is therefore reported separately
rather than counted as further calibration cases; the calibration remains 62
cases.

## Boundaries

- 62 cases, 35 of them on real reads. This is a calibration of what the tiers
  mean on these cases, not an estimate of exact-reconstruction sensitivity
  across taxa.
- The five closed genomes are one *H. pylori* collection sequenced on one
  platform. The depth series varies depth within those libraries; it does not
  vary read length, insert size, or error profile.
- `supported` spans few distinct loci, and no case where a single-copy locus is
  both close to the bait and genuinely hard to assemble. Acceptance within that
  class is therefore an upper bound.
- Truth for stages B to D is the RefSeq annotation of each closed genome. A
  boundary that differs from the biological one appears here as an inexact
  match, which is why the relation to the annotation is reported separately
  from the strict verdict: `relation_to_annotation` distinguishes disagreeing
  bases from an agreeing sequence that ends elsewhere.
- Stage A inherits the assembly and read simulation of the benchmarks it reads;
  it recomputes the calibration from their deposited verdicts rather than
  re-running them.
