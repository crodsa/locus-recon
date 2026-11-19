# Benchmark reproduction record

The compact files in this directory record the behaviour of the bundled synthetic
benchmark. They are reference results for this deterministic panel, not universal
estimates of biological performance.

The run recorded here was made from a clean results directory with Locus-Recon 1.0.1 on 2026-11-18.

```bash
python validation/run_mock_validation.py \
  --force \
  --threads 2 \
  --memory-gb 4
```

Environment:

- blastn: blastn: 2.16.0+
- bwa_mem2: 2.2.1
- python: Python 3.11.13
- samtools: samtools 1.22.1
- spades: SPAdes genome assembler v4.2.0
- package source digest (`source_sha256`, see `validation/common.py`):
  `2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0`

Outcome:

- 13 samples evaluated;
- workflow status: 1 `SKIP`, 12 `SUCCESS`;
- core criteria passed: 6 of 6 (`all_core_criteria_passed`: true);
- `baseline_intact_A`: status `SUCCESS`, tier `HIGH`, disposition `PASS`, identical to the truth true, remap depth 39.2x, mixed sites 0 (median alternative fraction 0.0), strand-biased sites 0, ambiguity flagged false;
- `fragmented_mid_A`: status `SUCCESS`, tier `HIGH`, disposition `PASS`, identical to the truth true, remap depth 38.7x, mixed sites 0 (median alternative fraction 0.0), strand-biased sites 0, ambiguity flagged false;
- `fragmented_terminal_A`: status `SUCCESS`, tier `HIGH`, disposition `PASS`, identical to the truth true, remap depth 38.6x, mixed sites 0 (median alternative fraction 0.0), strand-biased sites 0, ambiguity flagged false;
- `depth_03x_A`: status `SUCCESS`, tier `SUSPECT`, disposition `HOLD`, identical to the truth true, remap depth 3.2x, mixed sites 0 (median alternative fraction 0.0), strand-biased sites 0, ambiguity flagged false;
- `depth_05x_A`: status `SUCCESS`, tier `SUSPECT`, disposition `HOLD`, identical to the truth true, remap depth 4.8x, mixed sites 0 (median alternative fraction 0.0), strand-biased sites 0, ambiguity flagged false;
- `depth_10x_A`: status `SUCCESS`, tier `LOW`, disposition `REVIEW`, identical to the truth true, remap depth 9.7x, mixed sites 0 (median alternative fraction 0.0), strand-biased sites 0, ambiguity flagged false;
- `depth_20x_A`: status `SUCCESS`, tier `MEDIUM`, disposition `REVIEW`, identical to the truth true, remap depth 19.3x, mixed sites 0 (median alternative fraction 0.0), strand-biased sites 0, ambiguity flagged false;
- `mixture_05pct_B`: status `SUCCESS`, tier `SUSPECT`, disposition `HOLD`, identical to the truth true, remap depth 58.0x, mixed sites 2 (median alternative fraction 0.0588), strand-biased sites 0, ambiguity flagged false;
- `mixture_10pct_B`: status `SUCCESS`, tier `SUSPECT`, disposition `HOLD`, identical to the truth true, remap depth 58.0x, mixed sites 22 (median alternative fraction 0.0902), strand-biased sites 0, ambiguity flagged false;
- `mixture_20pct_B`: status `SUCCESS`, tier `SUSPECT`, disposition `HOLD`, identical to the truth true, remap depth 58.0x, mixed sites 25 (median alternative fraction 0.1875), strand-biased sites 0, ambiguity flagged false;
- `mixture_50pct_B`: status `SUCCESS`, tier `SUSPECT`, disposition `HOLD`, identical to the truth true, remap depth 57.8x, mixed sites 25 (median alternative fraction 0.4603), strand-biased sites 0, ambiguity flagged true;
- `paralog_control`: status `SUCCESS`, tier `SUSPECT`, disposition `HOLD`, identical to the truth false, remap depth 110.8x, mixed sites 14 (median alternative fraction 0.4718), strand-biased sites 0, ambiguity flagged true;
- `negative_control`: status `SKIP`, tier `none`, disposition `HOLD`, identical to the truth false, remap depth 0.0x, mixed sites 0 (median alternative fraction 0.0), strand-biased sites 0, ambiguity flagged false;

Timing: the 12 samples that entered read recruitment took a median of 5.20 s each (range 5.0 to 5.8 s), 63.9 s in total. These figures describe this run on this machine.

A second run of the same command was compared file by file with this one: 157 of 228 kept files were byte-identical (`evidence/mock_run_and_repeat_checksums.tsv` in the regeneration package lists every file).

In `locus_recon_report_mockLocus.tsv`, `allele_file` is shown relative to the results directory; the pipeline itself writes absolute paths.
