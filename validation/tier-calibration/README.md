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

62 cases with known truth, from four stages.

| case class | reached the top tier | top-tier precision | false accepts | withheld but exact |
|---|---|---|---|---|
| `supported` | **21 / 21** | **1.00** | 0 | 0 |
| `divergent` | 0 / 5 | — | 0 | 0 |
| `low-depth` | 0 / 10 | — | 0 | 9 |
| `truncated` | 0 / 11 | — | 0 | 11 |
| `multi-copy` | 0 / 5 | — | 0 | 4 |
| `degraded` | 0 / 10 | — | 0 | 6 |
| **all** | **21 / 62** | **1.00** | **0** | **30** |

Four statements follow, and they are the ones to quote.

1. **The top tier is reachable and it is precise.** Every case where the locus
   was single copy, within 95% of the bait, intact and adequately covered
   reached `HIGH`, and every result in `HIGH` matched truth exactly. There were
   no false accepts in 62 cases.
2. **An overall acceptance rate is not a performance measure.** It is 21 of 62
   here only because 41 of the cases were built or subsampled to be refused. On
   a dataset of fragmented or hypervariable targets the rate is expected to
   approach zero, and that is the designed behaviour rather than a failure.
3. **A withheld result is not a wrong result.** 30 of the 41 withheld cases
   reconstructed truth exactly over the span they reported. The tier states
   what the reads establish, not what the sequence happens to be, so a withheld
   result means "not established", never "incorrect".
4. **The one case with a wrong base was withheld.** Across all 62 cases exactly
   one reported sequence disagrees with truth in its bases: *gyrB* in Hpfe0002
   at 8x, at 99.96% over 2,322 bp, one substitution. It was reported at `LOW`
   with `REDUCED_REMAP_DEPTH (mean=5.5x)`, `PATCHY_REMAP_SUPPORT (34.3% <5x)`
   and `ELEVATED_UNCERTAIN_BASES (6.66% interior)`. Every other disagreement in
   the set is a boundary difference, not a base difference.

## Stage A: deposited benchmarks (27 cases, deterministic)

Recomputed from `../reference_results/validation_results.tsv` (13-sample mock
truth set) and `../completeness-benchmark/completeness_benchmark_results.tsv`
(14 constructed completeness cases), both of which carry a per-case truth
verdict alongside the tier. Per-case table:
`tier_calibration_benchmarks.tsv`.

Six of six `supported` cases reached `HIGH` with exact truth, including the two
where the locus was broken across contigs and rejoined. All eleven `truncated`
cases were reported at `MEDIUM` or `SUSPECT` and reconstructed exactly over the
recovered span. The paralogue control and the target-negative genome are the
only cases here that do not match truth, and both were withheld.

## Stage B: closed-genome audit on real reads (10 cases)

Five closed *Helicobacter pylori* chromosomes with their public Illumina runs.
The bait is a fixed external catalogue from strain 26695; sample-specific truth
is never used as a bait. Each draft is assembled with SPAdes `--isolate`, each
locus reconstructed, and each candidate compared to the annotated gene of its
own closed genome. Accessions: `closed_genome_manifest.tsv`; per-case results:
`tier_calibration_closed_genomes.tsv`.

All five *gyrB* loci reached `HIGH` and matched the annotated gene exactly
(2,322 bp) while the nearest catalogue allele sat at 96.2-96.7% identity, so
catalogue distance does not by itself suppress the tier. All five two-copy 23S
loci were withheld at `SUSPECT`; four are nevertheless exact, and Hpfe0006 is
truncated at 2,570 bp of 2,887.

## Stage C: divergent single-copy loci (15 cases)

The same drafts, with loci chosen by measurement rather than by hand. The
shortlist is prespecified in `CANDIDATE_LOCI`: single-copy *H. pylori* genes
with no known paralogue family, spanning housekeeping conservation to the
mosaic virulence genes. A locus is eligible only if the bait assembly and all
five closed genomes annotate it exactly once, and the three eligible loci
furthest from the bait allele are used. Every candidate and its per-genome copy
count is deposited in `locus_selection.tsv`; per-case results in
`tier_calibration_divergent_loci.tsv`.

| locus | mean identity to the bait | bait length | reached `HIGH` | truth |
|---|---|---|---|---|
| *ureB* | 97.2% | 1,710 bp | 5 / 5 | exact |
| *glmM* | 96.0% | 1,338 bp | 5 / 5 | exact |
| *cagA* | **88.7%** | 3,561 bp | 0 / 5 | identical over the overlap; boundary differs by -12 to +33 bp |

*cagA* is the informative case. At 88.7% identity to the one-allele bait the
tier falls to `MEDIUM` in three genomes and `SUSPECT` in two, and in all five
the reported sequence is **identical to the annotated allele over the whole
overlap** while ending somewhere else than the annotation does. So the
withholding here is not a response to wrong bases; it is a response to an
uncertain boundary in a locus whose 3' end carries a variable repeat, which is
the same regime as a hypervariable surface-antigen gene. Three candidates are
shorter than the annotation and two are longer.

Two loci in the shortlist were dropped for lack of a comparable annotation
rather than for biology: *vacA* is annotated under that gene name in the bait
assembly but in none of the five closed genomes, and *katA* and *flaA* appear
under that name in none of the annotations. The counts are in
`locus_selection.tsv`.

## Stage D: depth series (10 cases)

*gyrB*, the locus accepted in all five genomes at full depth, reconstructed
again from reads subsampled to ~15x and ~8x. Subsampling keeps every k-th pair,
which is deterministic and preserves the pairing, and the draft is reassembled
at each depth so the loss is felt by the assembly as well as by the
reconstruction. Per-case results: `tier_calibration_depth_series.tsv`.

| depth | tiers | truth |
|---|---|---|
| 576-923x (full) | `HIGH` x 5 | 5 exact |
| ~15x | `MEDIUM` x 3, `LOW` x 2 | 5 exact |
| ~8x | `LOW` x 4, `SUSPECT` x 1 | 4 exact, 1 with one wrong base |

The response is graded and monotone in every genome: no library that lost depth
kept its tier, and none was accepted. The reported sequence stays exact down to
8x in four of five genomes, which is the point of the exercise - the tool is
withholding sequence that happens to be right, because at 5-7x remapped depth
the reads do not establish it. The fifth genome is where that caution earns
itself: one substitution, withheld.

**The arms differ in assembler mode.** The full-depth drafts use SPAdes
`--isolate`, as the original audit did; the subsampled drafts use SPAdes'
default mode, because `--isolate` is documented for high-coverage isolate data.
Each arm therefore uses the settings one would actually use at that depth, and
the comparison across depths includes that difference.

## Command

```bash
# Stage A only: deterministic, no network or external tools
python validation/run_tier_calibration.py --outdir validation/tier-calibration

# All four stages
python validation/run_tier_calibration.py \
  --outdir validation/tier-calibration \
  --manifest validation/tier-calibration/closed_genome_manifest.tsv \
  --workdir /scratch/tier-calibration --threads 24 \
  --divergent-loci 3 --depth-series 15,8
```

## Boundaries

- 62 cases at six loci in three species, 35 of them on real reads. This is a
  calibration of what the tiers mean on these cases, not an estimate of
  exact-reconstruction sensitivity across taxa.
- The five closed genomes are one *H. pylori* collection sequenced on one
  platform. The depth series varies depth within those libraries; it does not
  vary read length, insert size, or error profile.
- `supported` contains 21 cases but only four distinct loci, and no case where a
  single-copy locus is both close to the bait and genuinely hard to assemble.
  Acceptance within that class is therefore an upper bound.
- Truth for stages B to D is the RefSeq annotation of each closed genome. A
  boundary that differs from the biological one appears here as an inexact
  match, which is why the relation to the annotation is reported separately
  from the strict verdict: `relation_to_annotation` distinguishes disagreeing
  bases from an agreeing sequence that ends elsewhere.
- Stage A inherits the assembly and read simulation of the benchmarks it reads;
  it recomputes the calibration from their deposited verdicts rather than
  re-running them.
