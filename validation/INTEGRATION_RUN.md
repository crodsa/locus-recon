# Integration run with the external tools

The 221 unit tests under `tests/` replace SPAdes, BWA-MEM2, samtools and BLAST+ with
stubs. They check the logic of the package and never start those programs. This file
documents one run of the complete workflow with the real programs, on the smallest
data set of the repository: the 13 simulated samples of `validation/mock_dataset/`.
The run can be repeated in about two minutes, and its outputs compared with the
checksums below.

## Environment of the recorded run

| Component | Version |
|---|---|
| locus-recon | locus-recon 1.0.1 |
| python | Python 3.11.13 |
| spades | SPAdes genome assembler v4.2.0 |
| bwa-mem2 | 2.2.1 |
| samtools | samtools 1.22.1 |
| blastn | blastn: 2.16.0+ |
| makeblastdb | makeblastdb: 2.16.0+ |
| biopython | 1.87 |
| numpy | 1.26.4 |
| platform | Linux 6.8.0-149-generic x86_64 |
| cpu | Intel(R) Core(TM) i7-14700 |

Package source digest (`source_sha256` of `validation/common.py`): `2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0`.
Git commit of the tree: `884cdf3c84f17dc9007576f26c6b431044dd35e1`.

## Commands

```bash
mamba env create -f environment.yml
mamba activate locus-recon
python -m pip install -e .
python validation/run_mock_validation.py \
    --results-dir work/mock_results --threads 2 --memory-gb 4 --force
```

The runner checks the data set, starts the workflow with the command below and then
evaluates the result against the simulated truth (`validation/evaluate_mock_results.py`):

```bash
python -m locus_recon.cli --samplesheet validation/mock_dataset/samples.tsv --main-output-dir work/mock_results --bait validation/mock_dataset/mock_locus_alleles.fasta --locus mockLocus --threads 2 --parallel-samples 1 --memory-per-sample 4 --min-bitscore-margin 80 --min-remap-depth 3 --min-remap-breadth 90 --min-base-quality 20 --min-mapping-quality 20 --mixture-min-fraction 0.05 --mixture-min-sites 2 --mixture-min-alt-depth 2 --no-progress
```

For every sample the workflow runs `makeblastdb` and `blastn` (discovery), `bwa-mem2 index` and
`bwa-mem2 mem` with `samtools` (mapping to the draft and recruitment of reads), `spades.py`
(local assembly), `blastn` (location and validation of the candidate) and `bwa-mem2` with
`samtools mpileup` (per-base support of the candidate).

## Expected result

`validation_summary.json`: `all_core_criteria_passed` is true; 6 of 6 criteria are met.

| Sample | Status | Tier | Disposition | Length (bp) | Identical to the truth | Mixed sites | SHA-256 of `<sample>_mockLocus_reconstructed.fasta` |
|---|---|---|---|---|---|---|---|
| `baseline_intact_A` | SUCCESS | HIGH | PASS | 507 | true | 0 | `647d2ade07b80da4154787c07cbd76971edf4a61d0d4dffd371f955bff669a0c` |
| `fragmented_mid_A` | SUCCESS | HIGH | PASS | 507 | true | 0 | `90b707ed72220e180cc7fe78175c3789f5fb51f5def8360ce180d926a69e3716` |
| `fragmented_terminal_A` | SUCCESS | HIGH | PASS | 507 | true | 0 | `83e2977da89b6366baace0c755fa722274726b0dbc8673eec67e68ba695fdaee` |
| `depth_03x_A` | SUCCESS | SUSPECT | HOLD | 507 | true | 0 | `ddf279a459febdb60a3139a83ceb4f356c7923cd8af396e886ba91fe7dc52eb4` |
| `depth_05x_A` | SUCCESS | SUSPECT | HOLD | 507 | true | 0 | `24909a1132c2a27eb5ad2eb0d0bf6a2bfcc6360d281793584e831e3eddef1f31` |
| `depth_10x_A` | SUCCESS | LOW | REVIEW | 507 | true | 0 | `dde4a7bcc3116307ae49ad5cf505b4131bfaeb3f63ebbb36955d552d69ca5b38` |
| `depth_20x_A` | SUCCESS | MEDIUM | REVIEW | 507 | true | 0 | `80fdd4c28d4079fd0d37d6527856361f5f26db5b52e355677ce716dcb998dc6d` |
| `mixture_05pct_B` | SUCCESS | SUSPECT | HOLD | 507 | true | 2 | `a4f67979b1af9d2ae477591c64555ab53656e89c701f84c601837676a65ef9cd` |
| `mixture_10pct_B` | SUCCESS | SUSPECT | HOLD | 507 | true | 22 | `159879d97d69d77f4eba8689e54ed7ed9cfc5088fd99eb9fca2b00c71c152f69` |
| `mixture_20pct_B` | SUCCESS | SUSPECT | HOLD | 507 | true | 25 | `9d11876bddb7adc089ec4ec32b5d172899d67e98ce685905a5c0500ffd94857a` |
| `mixture_50pct_B` | SUCCESS | SUSPECT | HOLD | 507 | true | 25 | `85a70da2ea654bb76b653619afe46eddc67f521139c4f7639842857fc019e642` |
| `paralog_control` | SUCCESS | SUSPECT | HOLD | 507 | false | 14 | `abb6c15a1527f5818c0ed23749a01004c460bc33a1a92825f738742105d05b40` |
| `negative_control` | SKIP | — | HOLD | 0 | false | — | no sequence |

Files of the results directory that summarise the run:

| File | Bytes | SHA-256 | Same bytes in a second run |
|---|---|---|---|
| `mockLocus_accepted_alleles.fasta` | 2109 | `be6a4bc860964622b0f7adf699b6bd52b02107c7c9fde605a6acaafbfdb974f4` | yes |
| `mockLocus_review_required_candidates.fasta` | 1378 | `e96ab116f8ca111bad6f0e5e77677afd5d5f87db4b78eaef667365d084a04add` | yes |
| `mockLocus_hold_candidates.fasta` | 4863 | `2844bc4c40f1837909601720129ece9db95a30ee08d779fc090f5acc4fb1ff6f` | yes |
| `mockLocus_reconstructed_candidates.fasta` | 8350 | `24b96160651d86f21c73f3e6952a8eb66b3261a9643aa6aa855c428eff1780f8` | yes |
| `validation/validation_summary.json` | 1174 | `c1654b26e0e03ca0b924fa5e2781c2600537076f18318792cbcaa2081dc9d666` | yes |
| `validation/validation_results.tsv` | 3096 | `c1d5865aeb718722dac6bf404832e6e9aa82eb4e9913633e54f2f130f66769d5` | no |
| `locus_recon_report_mockLocus.tsv` | 9439 | `858c37a23f07b3c9e484e5d75781479620f4136b1672c5638915d7f0928cb498` | no |

## What is reproducible to the byte

The command was run twice on the same machine (67 s and 66 s). Of the 228 files kept from the first run, 157 had the same bytes in the second.

| Kind of file | Same bytes | Different bytes |
|---|---|---|
| `allele_region.bed` | 11 | 1 |
| `allele_remap.flagstat.txt` | 12 | 0 |
| `allele_remap.support.tsv` | 12 | 0 |
| `allele_remap_bwa.log` | 0 | 12 |
| `assembly_graph_after_simplification.gfa` | 11 | 1 |
| `blast_step1.tsv` | 13 | 0 |
| `blast_step5.tsv` | 11 | 1 |
| `blast_validation.tsv` | 12 | 0 |
| `bwa_pipe.log` | 0 | 12 |
| `contigs.fasta` | 11 | 1 |
| `locus_recon_batch_mockLocus.log` | 0 | 1 |
| `locus_recon_report_mockLocus.tsv` | 0 | 1 |
| `mapping_flagstat.txt` | 12 | 0 |
| `mockLocus_accepted_alleles.fasta` | 1 | 0 |
| `mockLocus_hold_candidates.fasta` | 1 | 0 |
| `mockLocus_reconstructed_candidates.fasta` | 1 | 0 |
| `mockLocus_review_required_candidates.fasta` | 1 | 0 |
| `params.txt` | 0 | 12 |
| `qc_report.txt` | 12 | 0 |
| `recon.log` | 0 | 13 |
| `reconstructed.fasta` | 12 | 0 |
| `run_manifest.json` | 0 | 1 |
| `scaffolds.fasta` | 11 | 1 |
| `spades.log` | 0 | 12 |
| `target_regions.bed` | 12 | 0 |
| `validation_results.tsv` | 0 | 1 |
| `validation_summary.json` | 1 | 0 |
| `workflow_provenance.json` | 0 | 1 |

The candidate sequences, the QC reports, the per-base tables (`allele_remap.support.tsv`) and
`validation_summary.json` are the files to compare. Logs, the batch report, `validation_results.tsv`,
`run_manifest.json` and `workflow_provenance.json` carry dates, run times or absolute paths and differ
between runs. SPAdes writes its own paths and times in `params.txt` and `spades.log`.
In 1 sample(s) (`depth_10x_A`) the local SPAdes assembly itself differed between the two runs, and with it the files that describe the contig; the candidate sequence, its QC report and its per-base table were the same.

To compare a new run with the recorded one:

```bash
cd work/mock_results
sha256sum */*_mockLocus_reconstructed.fasta mockLocus_*.fasta validation/validation_summary.json
```

The complete list of checksums of the recorded run is `evidence/mock_run_and_repeat_checksums.tsv` in the
regeneration package (`validation/regenerated-1.0.1/`), and the files themselves are under `evidence/mock/`.

## Relation to the unit tests

| Test run | Last line of pytest |
|---|---|
| `release_as_committed.pytest.txt` | 1 failed, 220 passed in 21.76s |
| `release_with_all_records.pytest.txt` | 1 failed, 220 passed in 21.85s |
| `release_with_all_records_and_test_update.pytest.txt` | 221 passed in 25.53s |
| `release_with_liba6656_record.pytest.txt` | 221 passed in 23.33s |

All external programs are replaced with stubs in those tests. A change that breaks the call of an
external program, or that depends on its real output, is therefore visible only in a run such as this
one. The other analyses under `validation/` are also runs with the real programs; they need public
reads of 0.2 to 1 GB per sample and take hours.
