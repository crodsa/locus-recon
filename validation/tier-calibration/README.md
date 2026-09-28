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
`tier_calibration_tables.md` by the run itself, and `tier_calibration.pdf`
(also `.svg`) shows the tier by case class, the top-tier rate with its exact
interval, and the depth series.

### All four stages

| case class | reached the top tier | top-tier precision (95% CI) | false accepts | withheld but exact |
|---|---|---|---|---|
| `supported` | 21 / 21 | 1.00 (0.84-1.00) | 0 | 0 |
| `divergent` | 0 / 5 | — | 0 | 0 |
| `low-depth` | 0 / 10 | — | 0 | 9 |
| `truncated` | 0 / 11 | — | 0 | 11 |
| `multi-copy` | 4 / 5 | 1.00 (0.40-1.00) | 0 | 0 |
| `degraded` | 0 / 10 | — | 0 | 8 |
| **all** | 25 / 62 | 1.00 (0.86-1.00) | 0 | 28 |

The stages on real reads were run twice from the same inputs, and the 35
per-case results of stages B to D were identical between the runs.

Four statements answer the question the deposit asks, and one prespecified
expectation was not met.

1. **Is the top tier reachable, and is it precise?** All 21 supported cases
   reached `HIGH` and every one matched truth exactly. All 25 results reported
   at `HIGH` were exact (precision 1.00, exact 95% interval 0.86-1.00), and
   there were no false accepts in 62 cases.
2. **An overall acceptance rate is not a performance measure.** Most cases
   were built or subsampled to be refused, so the rate is a property of the
   case mix. On a dataset of fragmented or hypervariable targets it is expected
   to approach zero, and that is the designed behaviour rather than a failure.
3. **A withheld result is not a wrong result.** 28 of the 37 withheld cases
   reconstructed truth exactly over the span they reported. The tier states
   what the reads establish, not what the sequence happens to be, so a
   withheld result means "not established", never "incorrect".
4. **Where the reported bases disagree with truth, what tier did they get?**
   Two reported sequences disagree with truth in their bases, and both were
   withheld: the mock paralogue control (5 of 507 bases, at the positions that
   tell the paralogue apart; `SUSPECT`) and *gyrB* of Hpfe0002 at ~8x (99.96%
   over 2,322 bp, one substitution; `LOW`, with
   `REDUCED_REMAP_DEPTH (mean=5.6x)`, `PATCHY_REMAP_SUPPORT (33.2% <5x)` and
   `ELEVATED_UNCERTAIN_BASES (6.66% interior)`). Every other reported
   sequence that is not exact (Hpfe0006 23S and the five *cagA* alleles)
   differs from truth only at its boundary, as `relation_to_annotation` in the
   per-case tables records.

**The prespecified expectation for `multi-copy` was not met.** The class table
above expected multi-copy loci to be withheld even when the consensus is exact.
Four of the five two-copy 23S loci reached `HIGH` instead, each identical to
the annotated copy. The two copies are identical, so every read supports every
reported base, and the tier states that support: it is not a copy-number
statement. Checked as a coding sequence, the 23S gene would be withheld by frame
flags that have nothing to do with copy number; reconstructed with
`--noncoding-locus`, as a non-coding locus must be, it is not withheld. Copy
number is the question of the depth and graph modules, which on these genomes
return a depth excess over one copy at every 23S locus and two graph contexts
at four of the five, with a lower bound of two at the fifth, where the
traversal met a cycle (`../low-copy-23S/`).

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
| *gyrB* loci at `HIGH` that match the annotated gene exactly | 5 of 5 |
| Identity of the *gyrB* reconstructions to the nearest catalogue allele | 96.2-96.7% |
| Two-copy 23S loci: tiers, and consensus exact to the annotation | `HIGH` x 4, all exact; `SUSPECT` x 1 (Hpfe0006, below) |
| Reconstructions that were not exact, and how they were reported | Hpfe0006 23S: identical to the annotated copy over 2,570 bp but truncated by 317 bp at a contig end; `SUSPECT` and `HOLD`, with `LOCUS_SPLIT_ACROSS_CONTIGS` naming `NODE_2_length_784_cov_947.458274` as carrying all 317 missing bases |

Catalogue distance does not by itself suppress the tier: every *gyrB* allele is
exact at `HIGH` while its nearest catalogue allele sits 3-4 points away.

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
the reported sequence is identical to the annotated allele over the whole
overlap while ending somewhere else than the annotation does. The two
candidates that run past the annotated boundary also carry internal stop
codons in their best frame and are held at `SUSPECT`; the three that end short
of it carry only a length deviation from a one-allele profile and are reported
at `MEDIUM`. The withholding is a response to an uncertain boundary in a locus
whose 3' end carries a variable repeat, not to wrong bases.

A locus is counted as `divergent` only when its measured identity to the bait
is below 95%; a selected locus that sits as close to the bait as the
housekeeping locus is counted as `supported`. Where a divergent locus is
withheld, `relation_to_annotation` separates a response to wrong bases from a
response to an uncertain boundary.

Three shortlisted loci are shorter than the 1,200 bp minimum in the bait
(*recA*, *ftsZ*, *rocF*), and three have no comparable annotation: *katA* and
*flaA* appear under that gene name in no annotation, and *vacA* in the bait
assembly but in none of the five closed genomes. They are excluded for those
reasons, not for biology; the counts are in `locus_selection.tsv`.

## Stage D: depth series (10 cases)

*gyrB*, the locus accepted at full depth, reconstructed again from reads
subsampled to ~15x and ~8x. Subsampling keeps every k-th pair, which is
deterministic and preserves the pairing, and the draft is reassembled at each
depth so the loss is felt by the assembly as well as by the reconstruction.
Per-case results: `tier_calibration_depth_series.tsv`.

| depth | tiers | truth |
|---|---|---|
| full (stage B) | `HIGH` x 5 | 5 exact |
| ~15x | `MEDIUM` x 5 | 5 exact |
| ~8x | `LOW` x 4, `SUSPECT` x 1 | 4 exact, 1 with one wrong base |

The response is graded and monotone in every genome: no library that lost depth
kept its tier, and none was accepted. The reported sequence stays exact down to
~8x in four of five genomes; the tool withholds sequence that happens to be
right, because at 5-7x remapped depth the reads do not establish it. The fifth
genome, Hpfe0002, carries one substitution at ~8x and was withheld at `LOW`.
Hpfe0001 at ~8x is held at `SUSPECT` although exact, because one interior
position had no read support (`INTERNAL_ZERO_DEPTH_GAP`).

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
The reported tier is identical to the primary series in 10 of 10 cases, and the
agreement with truth is identical in 10 of 10, so the tier response follows
depth rather than the assembler setting.
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
