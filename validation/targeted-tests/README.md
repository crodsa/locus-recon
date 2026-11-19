# Targeted tests with simulated reads

This folder holds the targeted tests of section 3 of the pre-registered protocol
(`validation/external-evaluation/PROTOCOL.md`), run with Locus-Recon 1.0.1 (package digest
`2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0`). It is produced by
`validation/run_targeted_tests.py`. Every departure from the protocol is in `DEVIATIONS.md`.

The package was not modified. The runner calls the installed commands `locus-recon`,
`locus-recon-depth-ratio` and `locus-recon-copy-number` with default options (threads and memory
apart) and, for the block-size series only, the library function
`locus_recon.depth_ratio.estimate_locus_copy_number`.

## What was run

Every run is one simulated sample: one or two synthetic genomes with known alleles, paired
reads of 150 bp with substitution errors at a rate of 0.001, a draft assembly made with
`spades.py --isolate`, and the standard Locus-Recon workflow with a bait of one sequence. Tests 4
and 5 add the two copy-number commands. There are 790 runs; all 790 finished
(`targeted_manifest.tsv`, column `run_status`).

| Test | Design | Runs | Per-run table |
|---|---|---|---|
| 1 Mixtures with few sites | 1,500-bp locus; two alleles that differ at 1, 2, 5 or 10 sites; minor fraction 5, 10, 20, 35 or 50%; 60x; 10 replicates; 20 pure controls | 220 | `test1_mixtures_few_sites.tsv` |
| 2 Strand-bias alerts in pure samples | loci of 501, 1,500, 3,000 and 6,999 bp at 20, 60, 200 and 1,000x; 10 replicates | 160 | `test2_strand_bias_pure.tsv` |
| 3 Ends, indels and pseudogenes | 1,500-bp locus; nine sample alleles (control; a substitution within 5, 20 or 50 bp of each end; 3-bp deletion; 12-bp insertion; 1-bp deletion; two compensating indels 30 bp apart; premature stop codon); 10, 20 and 60x; 5 replicates | 135 | `test3_ends_indels_pseudogenes.tsv` |
| 4 Graph contexts | 3,000-bp locus; 17 scenarios with 1 to 4 copies; two backgrounds (`plain`, `decoy`); 60x; 5 replicates | 170 | `test4_graph_contexts.tsv` |
| 5 Overlapping spans | 3,000-bp locus; one copy, two identical copies, and two divergent copies each kept in part; discovery intervals that share 0, 10, 25, 45, 65, 83 or 90% of the query; 60x; 5 replicates | 105 | `test5_overlapping_spans.tsv` |
| Block size | bootstrap blocks of 100, 250, 500 and 1,000 bp for every run of tests 4 and 5 | 275 x 4 | `block_size_tests4_5.tsv` |

The truth of every run (alleles, variant positions, copy number, flank contexts, seeds) is in
`targeted_manifest.tsv` and in the per-run tables. The true alleles, the baits and the returned
sequences are in `sequences/`.

Vocabulary: accepted = `PASS`, sent to review = `REVIEW`, held = `HOLD`. "Exact" means that the
returned sequence is identical, base for base and over its whole length, to a true allele of
the sample.

## Results at a glance

Every number below is in `targeted_summary.tsv` (rows whose `test` ends in `-total` are marginal
totals) or in the table named beside it.

**Test 1.** A mixture that differs at one site was never held as a mixture (0 of 50) and was
always accepted (50 of 50); 41 of the 50 accepted sequences carried
`MIXED_SITE_BELOW_MIXTURE_CALL` and 9 carried no mixed-site flag (7 at 5% and 2 at 20%; in
the latter two the site is the third and the fifth base from the 3' end). With two sites, 2, 9, 10, 10 and 10 of 10
replicates were held as mixture at 5, 10, 20, 35 and 50%; with five sites 6, 10, 10, 10 and 10;
with ten sites 9, 10, 10, 10 and 10 (`test1_detection_limit.tsv`). At 5%, 6 of 10 two-site and
3 of 10 five-site mixtures were accepted. No accepted sequence was inexact (0 of 59 accepted
mixtures). The 24 returned sequences that match neither allele were all held; 23 of them are 50%
mixtures of two or more sites. The 20 pure controls were accepted, exact and without flag
(`test1_cells.tsv`).

**Test 2.** `STRAND_BIASED_ALT_SUPPORT` was raised in 4 of 160 pure samples: 3 of 10 at
6,999 bp and 20x, 1 of 10 at 501 bp and 1,000x, and none in the other 14 cells. The four
samples were held although their sequence was exact. One further pure sample (3,000 bp, 20x)
carried `MIXED_SITE_BELOW_MIXTURE_CALL`. All 160 sequences were exact. At 20x, 28 of 40
samples were sent to review, all with `REMAP_DEPTH_BELOW_HIGH`, and 9 were accepted; at 60x
and 200x all 80 were accepted (`test2_cells.tsv`).

**Test 3.** 133 of 135 sequences were exact, and all 45 at 60x. The two inexact sequences (one
wrong base each, both at 10x) were held. No inexact sequence was accepted. Substitutions near
the ends were recovered in 44 of 45 runs. The 12-bp insertion was never accepted (0 of 15,
all exact): its length is flagged against a bait profile of one allele (`LENGTH_ANOMALOUS`).
The frame flags agreed with the truth in 134 of 135 runs. All 37 samples whose true reading
frame is disrupted (15 with the 1-bp deletion, 15 with the premature stop codon and 7 of the
compensating indels) were held; 36 of these sequences were exact, so a correct pseudogene
allele is not accepted. Uncertain bases were at the ends: the median number in the interior
was 0 in 26 of the 27 cells and 1 in the other, and the median over the whole sequence ranged
from 2 to 74 (`test3_cells.tsv`).

**Test 4.** With identical copies that the assembler collapsed (decoy background, 2 to 4
dispersed copies, shared flanks of 500 and 2,000 bp, two copies in a strain mixture; 30 runs)
the graph count equalled the true copy number, depth called more than one copy and the
consensus was the true integer in 30 of 30; the 10 runs with 3 and 4 identical copies in the
plain background gave the same. In single-copy samples the count was 1 in 20 of 20
controls and flank-repeat samples; with a 50:50 strain mixture that differs at one flank site
the two sides gave 2 and 1 contexts and no count was returned (10 of 10,
`ASYMMETRIC_CONTEXTS`). Unfavourable results: (i) when the draft held two complete identical
copies (plain background, 20 runs) the graph count was 2 but the depth ratio was 0.000, the
depth call was `SINGLE_COPY_COMPATIBLE` and the consensus was `INDETERMINATE` in 20 of 20;
(ii) a tandem duplicate of identical copies was not quantified in any of 10 runs (graph lower
bound 1, dosage 0.000, consensus `INDETERMINATE`); (iii) with copies 1% divergent the graph gave
no count in 69 of 70 runs (`INSUFFICIENT_TARGET_COVERAGE`) and the depth call was
`MULTICOPY_DEPTH` in only 28 of 70, although the dosage was within 0.5 of the truth in 66 of 70;
the consensus had the status `REVIEW` in 66 of 70. No consensus with the status `SUPPORTED` was wrong (0 of 170 with
the graph, 0 of 170 without) (`test4_cells.tsv`).

**Test 5.** For drafts in which the overlapped bases are present twice with identical sequence,
the dosage separated one copy from two at 0, 10 and 25% of overlap (medians 1.03, 1.05 and 1.10
against 1.93, 2.15 and 2.27) and was marked `REVIEW_OVERLAPPING_SPANS` whenever spans
overlapped. At 45% both were `LOWER_BOUND` (0.38 against 0.76), and at 65, 83 and 90% both
returned 0.000 with `LOWER_BOUND` and an `INDETERMINATE` consensus: at these overlaps the
dosage does not separate one copy from two. The bootstrap intervals of the two copy numbers
were disjoint only at 0 and 10% (`test5_separation.tsv`). At 25% the depth call for two copies
was `SINGLE_COPY_COMPATIBLE` in 5 of 5 and the consensus with the graph was
`EVIDENCE_CONFLICT`. In the added arm with two divergent copies, each kept in part, the median
dosage was 1.958, 2.136, 2.240, 1.800, 1.705, 1.769 and 1.853 at 0, 10, 25, 45, 65, 83 and 90%,
always with `REVIEW_OVERLAPPING_SPANS` when spans overlapped; at 83% it was 1.769 against a
`geometry_expected_dosage` of 1.83, and the depth call was `SINGLE_COPY_COMPATIBLE` from 25%
upwards. Without the graph the consensus was never `SUPPORTED` when spans overlapped
(`test5_cells.tsv`).

**Block size.** `block_size_summary.tsv` gives, for each test, background and block size, the
interval width and how often the dosage interval contains the true copy number;
`block_size_cells.tsv` gives the same by scenario. In test 4 with the decoy background the
interval contained the truth in 52, 57, 55 and 49 of 85 runs with blocks of 100, 250, 500 and
1,000 bp, and the depth call differed from the 500-bp call in 3, 2, 0 and 5 runs.

## Files

| File | Content |
|---|---|
| `README.md`, `DEVIATIONS.md` | this description; every departure from the protocol |
| `PROVENANCE.json` | versions, package digest, command templates, seeds, simulation parameters |
| `targeted_manifest.tsv` | the 790 runs: design factors, truth, seeds, status, directory on the host |
| `targeted_summary.tsv` | headline numbers in long format (`test`, `cell`, `metric`, `value`, `denominator`, `source_file`) |
| `test1_mixtures_few_sites.tsv`, `test1_cells.tsv`, `test1_detection_limit.tsv` | test 1: every run; counts per cell; detection limit by number of sites |
| `test2_strand_bias_pure.tsv`, `test2_cells.tsv` | test 2: every run; counts per length and depth |
| `test3_ends_indels_pseudogenes.tsv`, `test3_cells.tsv` | test 3: every run; counts per scenario and depth |
| `test4_graph_contexts.tsv`, `test4_cells.tsv` | test 4: every run; counts per background and scenario |
| `test5_overlapping_spans.tsv`, `test5_cells.tsv`, `test5_separation.tsv` | test 5: every run; counts per arm and overlap; separation of one copy from two |
| `block_size_tests4_5.tsv`, `block_size_cells.tsv`, `block_size_summary.tsv` | block-size series: every run and block size; per scenario; per test and background |
| `input_checksums.tsv` | SHA-256 and size of the reads, bait, truth, draft, graph and alignment of every run |
| `sequences/testN_truth_alleles.fasta`, `testN_baits.fasta`, `testN_returned.fasta` | true alleles (headers `run_id|allele`), baits and returned sequences, per test |
| `logs/` | log and timing of the full run, log of the four repeated runs, and the diff between the first and the delivered version of the runner |
| `pilot/` | tables of the three pilots (replicate 91; not part of the results) |

Columns of the per-run tables: design factors and truth first; `lr_*` are columns of the
Locus-Recon report (`lr_flag_names` lists the flag names, `lr_qc_flags` the full flag text);
`exact`, `matches_allele` and `cmp_*` compare the returned sequence with the truth; `dr_*` are
the fields of `locus-recon-depth-ratio`; `left_context_count` to `graph_context_flags` and
`cn_*` come from `locus-recon-copy-number` with the graph, `cn_nograph_*` from the same command
without the graph; `graph_vs_truth`, `cn_vs_truth` and `depth_call_agrees_with_truth` score
them against the true copy number.

## How to rerun

```bash
# environment.yml of the repository, with Locus-Recon 1.0.1 installed, plus numpy and biopython
python validation/run_targeted_tests.py \
    --work-dir /scratch/targeted_work \
    --out-dir validation/targeted-tests \
    --jobs 8 --threads-per-run 1 --memory-gb 8
```

The runner refuses to start if the installed package is not 1.0.1 with the frozen digest. All
seeds derive from the base seed 20261118 and the run identifier, so the same command simulates
the same samples; `input_checksums.tsv` allows the comparison. A run that has a `result.json`
in the work directory is not repeated. `--tests 1,3` restricts the tests, `--only RUN_ID,...`
restricts the runs and `--aggregate-only` rebuilds the tables from the stored records. The full
run took 1 h 1 min with eight single-thread runs in parallel (`logs/full_run_time.txt`); the
size of the work directory is in `logs/disk_usage.txt`.

## Limits

The reads carry only substitution errors, at the same rate on both strands and with constant
base quality, and the genomes are 100 kb of random sequence. The tests show how the decision
rules respond to sampling noise and to defined structures. They do not reproduce
platform-specific errors, G+C bias or the repeat content of a bacterial genome. Each cell has
5 or 10 replicates, so a count of 0 in a cell does not exclude a rate of several per cent.
