# Benchmark reproduction record

The compact files in this directory record the expected behaviour of the bundled
synthetic benchmark. They are reference results for this deterministic panel,
not universal estimates of biological performance.

The run recorded here was made from a clean results directory with
Locus-Recon 1.0.0 on 2026-09-25.

```bash
python validation/run_mock_validation.py \
  --force \
  --threads 2 \
  --memory-gb 4
```

Environment:

- Linux x86_64
- Python 3.11.13
- BWA-MEM2 2.2.1
- samtools 1.22.1
- SPAdes 4.2.0
- BLAST+ 2.16.0+
- package source digest (`source_sha256`, see `validation/common.py`):
  `3863abc32fb5446a72919c790cb2558217871410f441db1ecf17dde3de55cdcf`

Outcome:

- 13 samples evaluated;
- 12 `SUCCESS`, one `SKIP` (the target-negative control), zero execution
  failures;
- all six prespecified core criteria passed;
- exact truth for the intact and both fragmented controls, all three `HIGH`;
- exact truth at every depth of the sensitivity series, at observed remap depths
  of 3.2x, 4.8x, 9.7x and 19.3x for nominal 3x, 5x, 10x and 20x, with tiers
  `SUSPECT`, `SUSPECT`, `LOW` and `MEDIUM`: the sequence is right and the tier
  reports how thin its read support is;
- 5%, 10%, 20% and 50% mixtures detected, with 2, 22, 25 and 25 bidirectional
  mixed sites and median alternative fractions of 0.0588, 0.0902, 0.1875 and
  0.4603; the 5% mixture sits exactly at the two-site minimum, which is why the
  core criterion covers 10% and above;
- no strand-biased candidate site in any mixture sample;
- no mixture call in any successful pure control;
- close paralogue marked ambiguous;
- target-negative control not reconstructed.

Per-base support counts every read that passes the base- and mapping-quality
filters, including pairs the aligner could not mark as properly paired against a
507-bp reference (`samtools mpileup -A`), so observed remap depth tracks nominal
depth and both strands contribute at every site.

Timing, with two threads per sample, serial sample processing and a 4 GB SPAdes
memory cap: the 12 samples that entered read recruitment took a median of 5.65 s
each (range 5.4 to 6.7 s), 70.3 s in total; the target-negative sample ended at
discovery in 0.1 s. These figures reproduce this run on comparable hardware;
they are not a performance estimate for other loci, datasets or machines.

The generator seeds every sample and hashes every generated file, so the run is
re-executable and byte-comparable. Two independent runs of the same code
returned identical `validation_summary.json` files and identical per-sample
biological and QC metrics; only run dates and elapsed times differed. In the
committed `locus_recon_report_mockLocus.tsv`, `allele_file` is shown relative
to the results directory; the pipeline itself writes absolute paths.

Full local outputs are created under `validation/mock_results/` and are
intentionally excluded from source distributions because they include bulky BAM
and SPAdes intermediates. The deterministic generator and runner recreate them
from the committed truth set.
