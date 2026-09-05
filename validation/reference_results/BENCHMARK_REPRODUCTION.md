# Benchmark reproduction record

The compact files in this directory record the expected behaviour of the bundled
synthetic benchmark. They are reference results for this deterministic panel,
not universal estimates of biological performance.

The run recorded here was made from a clean results directory with the final
Locus-Recon 1.0 release on 2026-08-24.

```bash
python validation/run_mock_validation.py \
  --force \
  --threads 2 \
  --memory-gb 4
```

Environment:

- macOS arm64
- Python 3.14.0
- BWA-MEM2 2.2.1
- samtools 1.22.1
- SPAdes 4.2.0
- BLAST+ 2.16.0+

Outcome:

- 13 samples evaluated;
- 10 `SUCCESS`, three `SKIP`, zero execution failures;
- all six prespecified core criteria passed;
- exact truth for the intact and both fragmented controls;
- 3x and 5x depth controls skipped;
- exact truth at nominal 10x and 20x, at observed remap depths of 7.8x and
  14.4x;
- 5%, 10%, 20%, and 50% mixtures detected, with median alternative fractions of
  0.0526, 0.0851, 0.1984, and 0.4682 respectively;
- no mixture call in any successful pure control;
- close paralogue marked ambiguous;
- target-negative control not reconstructed.

Timing, with two threads per sample, serial sample processing and a 4 GB SPAdes
memory cap: the 12 samples that entered read recruitment took a median of 5.8 s
each (range 5.6 to 6.1 s), 69.8 s in total; the target-negative sample ended at
discovery in under 0.1 s. These figures reproduce this run on comparable
hardware; they are not a performance estimate for other loci, datasets or
machines.

The generator seeds every sample and hashes every generated file, so the run is
re-executable and byte-comparable. Reproducing it on the same software stack
returns a `validation_summary.json` identical to the committed one; per-sample
biological and QC metrics are identical too, and only run dates and elapsed
times differ.

Full local outputs are created under `validation/mock_results/` and are
intentionally excluded from source distributions because they include bulky BAM
and SPAdes intermediates. The deterministic generator and runner recreate them
from the committed truth set.
