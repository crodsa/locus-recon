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
| `supported` | single-copy locus, intact in the draft, adequate depth, pure culture | accept |
| `truncated` | the locus runs off a contig end by construction | tier falls with the size of the loss |
| `multi-copy` | the closed genome annotates more than one copy, so near-identical copies collapse in a short-read draft | withhold, even when the reported consensus is exact |
| `degraded` | depth below the stated floor, a controlled allele mixture, a paralogue, or a target-negative genome | withhold |

## Result

37 cases with known truth, from two independent stages.

| case class | reached the top tier | top-tier precision | false accepts | withheld but exact |
|---|---|---|---|---|
| `supported` | **11 / 11** | **1.00** | 0 | 0 |
| `truncated` | 0 / 11 | — | 0 | 11 |
| `multi-copy` | 0 / 5 | — | 0 | 4 |
| `degraded` | 0 / 10 | — | 0 | 6 |
| **all** | **11 / 37** | **1.00** | **0** | **21** |

Three statements follow, and they are the ones to quote.

1. **The top tier is reachable and it is precise.** Every case where the locus
   was single copy, intact and adequately covered reached `HIGH`, and every
   result in `HIGH` matched truth exactly. There were no false accepts.
2. **An overall acceptance rate is not a performance measure.** It is
   11 of 37 here only because 26 of the 37 cases were built to be refused. On a
   dataset of fragmented or hypervariable targets the rate is expected to be
   near zero, and that is the designed behaviour rather than a failure.
3. **A withheld result is not a wrong result.** 21 of the 26 withheld cases
   reconstructed truth exactly over the span they reported. The tier states
   what the reads establish, not what the sequence happens to be, so `HOLD`
   means "not established", never "incorrect".

## Stage A: deposited benchmarks (27 cases, deterministic)

Recomputed from `../reference_results/validation_results.tsv` (13-sample mock
truth set) and `../completeness-benchmark/completeness_benchmark_results.tsv`
(14 constructed completeness cases), both of which carry a per-case truth
verdict alongside the tier. Per-case table:
`tier_calibration_benchmarks.tsv`.

Six of six `supported` cases reached `HIGH` with exact truth, including the two
where the locus was broken across contigs and rejoined. All eleven `truncated`
cases were reported at `MEDIUM` or `SUSPECT` and reconstructed exactly over the
recovered span. The paralogue control and the target-negative genome were the
only cases that do not match truth, and both were withheld.

## Stage B: closed-genome audit on real reads (10 cases)

Five closed *Helicobacter pylori* chromosomes with their public Illumina runs.
The bait is a fixed external catalogue from strain 26695
(`bait_gyrB_26695.fasta`, `bait_23S_26695.fasta`); sample-specific truth is
never used as a bait. Each draft is assembled with SPAdes `--isolate`, each
locus reconstructed, and each candidate compared to the annotated gene of its
own closed genome. Accessions are in `closed_genome_manifest.tsv`; per-case
results in `tier_calibration_closed_genomes.tsv`.

| strain | locus | tier | matches truth | evidence |
|---|---|---|---|---|
| Hpfe0001 | *gyrB* | `HIGH` | yes | 2,322 bp exact, catalogue identity 96.4% |
| Hpfe0002 | *gyrB* | `HIGH` | yes | 2,322 bp exact, catalogue identity 96.2% |
| Hpfe0003 | *gyrB* | `HIGH` | yes | 2,322 bp exact, catalogue identity 96.6% |
| Hpfe0004 | *gyrB* | `HIGH` | yes | 2,322 bp exact, catalogue identity 96.3% |
| Hpfe0006 | *gyrB* | `HIGH` | yes | 2,322 bp exact, catalogue identity 96.7% |
| Hpfe0001 | 23S | `SUSPECT` | yes | 2,889 bp exact, two identical copies |
| Hpfe0002 | 23S | `SUSPECT` | yes | 2,887 bp exact, two identical copies |
| Hpfe0003 | 23S | `SUSPECT` | yes | 2,889 bp exact, two identical copies |
| Hpfe0004 | 23S | `SUSPECT` | yes | 2,889 bp exact, two identical copies |
| Hpfe0006 | 23S | `SUSPECT` | no | 2,570 bp of 2,887, truncated |

The five *gyrB* loci are the important row: the top tier was reached on real
reads, with truth known, while the nearest catalogue allele sat at
96.2-96.7% identity. Catalogue distance therefore does not suppress the tier,
which is the behaviour the confidence layer is designed to have: the tier
describes read support for the reported bases, not similarity to a reference.

The five 23S loci are the counterpart. Each is annotated twice in its own closed
genome, so a short-read draft collapses the copies and the reads cannot
establish which copy the consensus belongs to. The tool withholds all five, and
four of them are nevertheless exact. The fifth is truncated at 2,570 bp of
2,887 and is the single case in this deposit whose sequence does not match
truth; it was withheld.

## Command

```bash
# Stage A only: deterministic, no network or external tools
python validation/run_tier_calibration.py --outdir validation/tier-calibration

# Both stages
python validation/run_tier_calibration.py \
  --outdir validation/tier-calibration \
  --manifest validation/tier-calibration/closed_genome_manifest.tsv \
  --workdir /scratch/tier-calibration --threads 24
```

## Boundaries

- 37 cases at four loci in three species, 10 of them on real reads. This is a
  calibration of what the tiers mean on these cases, not an estimate of
  exact-reconstruction sensitivity across taxa.
- The `supported` class is small (11 cases) and contains no case where the
  locus was both single copy and genuinely difficult, so the acceptance rate
  within that class is an upper bound.
- Stage A inherits the assembly and read simulation of the benchmarks it reads;
  it recomputes the calibration from their deposited verdicts rather than
  re-running them.
- Truth for stage B is the RefSeq annotation of each closed genome. An
  annotation boundary that differs from the biological one would appear here as
  an inexact match rather than as an annotation difference.
