# External evaluation of Locus-Recon 1.0.1: results

This folder holds the results of section 2 of `validation/external-evaluation/PROTOCOL.md`, which was written
on 18 November 2026 before any of these analyses was run. Copies of the protocol, of the isolate table and of
the selection log are here as `protocol_PROTOCOL.md`, `protocol_isolates.tsv` and
`protocol_isolate_selection_log.tsv`. Every departure from the protocol is in `DEVIATIONS.md`. Read that file
with the numbers: one isolate has reads that do not match its deposited genome, 12 units have a bait that the
protocol rule does not find, and 3 units have no truth sequence.

The runner is `validation/run_external_evaluation.py`. It was run on 18 and 19 November 2026 on a Linux
workstation (28 cores) with at most 12 threads.

## What was run

- **Software under test:** Locus-Recon 1.0.1. The digest of `locus_recon/*.py` was checked at the start of
  every stage that uses the package, for the installed copy and for the copy in the repository, and both equal
  `2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0` (`PROVENANCE.json`, `frozen_package`). No file of the
  package and no threshold was changed.
- **Isolates and reads:** the 16 FDA-ARGOS isolates of `protocol_isolates.tsv`. Each Illumina run was downloaded
  from the European Nucleotide Archive, checked against the MD5 that the archive publishes, and subsampled to
  100x of the genome length with `seqtk sample -s 20261118` (`input_reads.tsv`). All 16 runs held more than
  100x. A second subsampling reproduced every file byte for byte.
- **Drafts:** SPAdes 4.2.0 `--isolate`, 8 threads, `contigs.fasta` (`isolates_and_drafts.tsv`).
- **Baits and truth:** one bait per locus from the RefSeq reference genome of the species (`baits.tsv`,
  `baits.fasta`); truth copies from the complete genome of the isolate by BLAST+ (`truth_copies.tsv`,
  `truth_sequences.fasta`).
- **Locus-Recon:** one run per unit with default options, `--threads 4`, and `--noncoding-locus` for the 23S and
  16S rRNA genes; `--cleanup` was not passed. Then `locus-recon-depth-ratio` and `locus-recon-copy-number` with
  the intervals of `blast_step1.tsv` and the graph of the draft.
- **Comparators:** ARIBA 2.14.7 and KMA 1.7.1 with default options on the same reads and the same baits, one
  catalogue of eight baits per isolate, through the functions of `validation/run_comparator_benchmark.py`.
  ARIBA accepted every bait as declared; no bait had to be declared non-coding.
- **Scoring:** `compare()` of `run_comparator_benchmark.py` for every returned sequence, and the decision rules
  of protocol 2.5. The ARIBA flag bits were confirmed with `ariba flag` on the host (`logs/host_probe_2026-11-18.txt`).

## Results in brief

Every number below is in the table named beside it. Intervals are 95% percentile intervals from 10,000
resamples of the 16 isolates (seed 20261118); no interval treats the units as independent.

### Units

128 units of isolate and locus were run. The bait rule of the protocol finds 116 of them (the protocol
set): 81 with one copy in the complete genome, 6 with several identical copies, 26 with several
different copies, and 3 without a truth copy because the bait aligns to the genome only in parts. The other
12 units have a bait found by a fallback rule and are reported apart as the supplementary set
(`external_units.tsv`, `baits.tsv`). No unit was a true negative control.

### One isolate has reads that do not match its genome

The draft of *Neisseria meningitidis* FDAARGOS_209 (run SRR4272082) aligns to 93.77% of the complete genome
GCF_002073575.2 at 97.35% identity with 45,838 substitutions, and only 67.84% of the draft aligns
(2,940,763 bp in 3,193 contigs for a genome of 2,181,232 bp). The other 15 drafts align to at least
98.85% of their genome at 99.97% identity or more, with at most 246 substitutions
(`isolates_and_drafts.tsv`; MUMmer dnadiff, not part of the protocol). The truth of this isolate therefore does
not describe its reads. All 16 isolates are kept in every table. The summaries are given for the protocol set
and, as an analysis added afterwards, for the protocol set without this isolate
("protocol, isolates whose draft agrees with the complete genome").

### Sequences against the truth, protocol set (113 units with a truth copy)

Identical to a truth copy over the whole length: Locus-Recon 98 of 113 (86.7%; interval over isolates 71.3 to 96.8%); KMA 99 of 113 (87.6%; interval over isolates 71.9 to 98.2%);
ARIBA 94 of 113 (83.2%; interval over isolates 68.4 to 93.2%). Without FDAARGOS_209: Locus-Recon 98 of 105, KMA 99 of 105, ARIBA 94 of 105
(`external_summary.tsv`).

Decision against result (`decision_tables.tsv`):

| Decision rule | Units | Accepted, exact | Accepted, not exact | Not accepted, exact | Not accepted, not exact |
|---|---|---|---|---|---|
| all strata with a copy: Locus-Recon, PASS | 113 | 89 | 6 | 9 | 9 |
| all strata with a copy: KMA, primary rule | 113 | 90 | 8 | 9 | 6 |
| all strata with a copy: KMA, secondary rule | 113 | 94 | 9 | 5 | 5 |
| all strata with a copy: ARIBA, protocol rule | 113 | 82 | 8 | 12 | 11 |
| one copy: Locus-Recon, PASS | 81 | 75 | 5 | 1 | 0 |
| one copy: KMA, primary rule | 81 | 73 | 5 | 2 | 1 |
| one copy: KMA, secondary rule | 81 | 74 | 5 | 1 | 1 |
| one copy: ARIBA, protocol rule | 81 | 74 | 6 | 0 | 1 |
| several identical copies: Locus-Recon, PASS | 6 | 5 | 1 | 0 | 0 |
| several identical copies: KMA, primary rule | 6 | 4 | 1 | 1 | 0 |
| several identical copies: KMA, secondary rule | 6 | 4 | 1 | 1 | 0 |
| several identical copies: ARIBA, protocol rule | 6 | 4 | 1 | 0 | 1 |
| several different copies: Locus-Recon, PASS | 26 | 9 | 0 | 8 | 9 |
| several different copies: KMA, primary rule | 26 | 13 | 2 | 6 | 5 |
| several different copies: KMA, secondary rule | 26 | 16 | 3 | 3 | 4 |
| several different copies: ARIBA, protocol rule | 26 | 4 | 1 | 12 | 9 |

The same table without FDAARGOS_209:

| Decision rule | Units | Accepted, exact | Accepted, not exact | Not accepted, exact | Not accepted, not exact |
|---|---|---|---|---|---|
| all strata with a copy: Locus-Recon, PASS | 105 | 89 | 0 | 9 | 7 |
| all strata with a copy: KMA, primary rule | 105 | 90 | 1 | 9 | 5 |
| all strata with a copy: KMA, secondary rule | 105 | 94 | 2 | 5 | 4 |
| all strata with a copy: ARIBA, protocol rule | 105 | 82 | 1 | 12 | 10 |
| one copy: Locus-Recon, PASS | 76 | 75 | 0 | 1 | 0 |
| one copy: KMA, primary rule | 76 | 73 | 0 | 2 | 1 |
| one copy: KMA, secondary rule | 76 | 74 | 0 | 1 | 1 |
| one copy: ARIBA, protocol rule | 76 | 74 | 1 | 0 | 1 |
| several identical copies: Locus-Recon, PASS | 5 | 5 | 0 | 0 | 0 |
| several identical copies: KMA, primary rule | 5 | 4 | 0 | 1 | 0 |
| several identical copies: KMA, secondary rule | 5 | 4 | 0 | 1 | 0 |
| several identical copies: ARIBA, protocol rule | 5 | 4 | 0 | 0 | 1 |
| several different copies: Locus-Recon, PASS | 24 | 9 | 0 | 8 | 7 |
| several different copies: KMA, primary rule | 24 | 13 | 1 | 6 | 4 |
| several different copies: KMA, secondary rule | 24 | 16 | 2 | 3 | 3 |
| several different copies: ARIBA, protocol rule | 24 | 4 | 0 | 12 | 8 |

- **Accepted and not identical (protocol 2.4, item 1):** Locus-Recon accepted 6 sequences that differ from
  the truth (`accepted_not_identical.tsv`). All 6 belong to FDAARGOS_209 (*gyrB*, *rpoB*, *recA*, *dnaK*, *atpD* and the
  16S rRNA gene; 3 to 86 substitutions, no gap). KMA returned the same six sequences and ARIBA the same five and
  a sixth that is 2 bp shorter, and both accepted them under their rules. In the other 15 isolates no sequence
  accepted by Locus-Recon differs from the truth (89 of 89 accepted are identical).
- **Proportion accepted (item 2):** Locus-Recon accepted 95 of 113 (84.1%; interval over isolates 76.3 to 90.8%): 80 of 81 one-copy
  units, 6 of 6 with several identical copies and 9 of 26 with several different copies.
- **Not accepted (item 3):** of the 18 units that Locus-Recon did not accept, 9 were complete and exact,
  3 incomplete and identical over their length, 6 wrong, and 0 had no sequence. Seventeen of the 18 are rRNA
  genes or *tuf* with several different copies; the eighteenth is the one-copy 16S rRNA gene of
  *M. tuberculosis* FDAARGOS_750, which is exact, was held at tier SUSPECT, and carries one flag,
  `STRAND_BIASED_ALT_SUPPORT` (1 of 1 candidate sites).
- **What the decision adds, and what it does not:** in units with several different copies Locus-Recon accepted
  9 sequences and all are exact, and it did not accept 8 exact ones. The primary rule of KMA accepted 15
  (2 not exact), its secondary rule 19 (3 not exact), and the rule of ARIBA 5 (1 not exact). Without
  FDAARGOS_209 the inexact accepted sequences in the whole set are 0 for Locus-Recon, 1 and 2 for the two rules of
  KMA and 1 for ARIBA, while Locus-Recon accepts 89 units, the rules of KMA 91 and 96, and ARIBA 83.
  In one-copy units the three tools behave alike. The decision of Locus-Recon did not detect the reads of
  FDAARGOS_209; neither did the rules of the other two tools.
- **Flags:** `STRAND_BIASED_ALT_SUPPORT` was raised in 1 of 81 one-copy units of the protocol set
  (1.2%). 7 of the 95 accepted sequences carry `MIXED_SITE_BELOW_MIXTURE_CALL`; they hold 7 mixed sites, and
  the truth copies differ at 5 of them (`mixed_sites.tsv`). The other two are a +TT insertion at 0.14 of the
  reads in the 23S rRNA gene of FDAARGOS_197 and a substitution at 0.06 of 36 reads in *recA* of FDAARGOS_1365.
- **Units without a truth copy:** the 23S rRNA gene of *H. influenzae* FDAARGOS_199 and of *S. enterica*
  FDAARGOS_54, and *recA* of *B. cereus* FDAARGOS_780. Every tool returned a sequence for each and none was
  accepted under any rule. Exactness is not defined for them.

### Copy number, protocol set

- **Threshold of 1.5:** 0 of 81 one-copy units were called `MULTICOPY_DEPTH`. Of the units with two or more
  copies, 21 of 32 (65.6%; interval over isolates 41.2 to 89.7%) were called `MULTICOPY_DEPTH`. The 11 that were not: FDAARGOS_57 16S (truth 6, 3 intervals, ratio 4.46); FDAARGOS_780 23S (truth 13, 3 intervals, ratio 9.33); FDAARGOS_780 16S (truth 13, 3 intervals, ratio 11.44); FDAARGOS_51 23S (truth 10, 6 intervals, ratio 8.13); FDAARGOS_51 16S (truth 10, 4 intervals, ratio 8.51); FDAARGOS_51 tuf (truth 2, 2 intervals, ratio 0.76); FDAARGOS_209 23S (truth 4, 6 intervals, ratio 1.83); FDAARGOS_209 16S (truth 4, 2 intervals, ratio 1.60); FDAARGOS_209 tuf (truth 2, 3 intervals, ratio 1.70); FDAARGOS_501 23S (truth 4, 6 intervals, ratio 0.28); FDAARGOS_501 16S (truth 4, 6 intervals, ratio 0.16).
  In all 11 the locus lies on two or more intervals of the draft; three belong to FDAARGOS_209.
- **Error of the ratio:** for one-copy units the median ratio is 1.087 (range 0.787 to 1.589) and the median
  absolute error 0.115. The ratio rises with G+C: the median relative error is +0.020 in genomes below
  40% G+C (37 units), +0.040 from 40 to 55% (26 units) and +0.205 at 55% or more (18 units)
  (`copy_number_error_by_factor.tsv`). Errors by copy number are in `external_summary.tsv`.
- **Bootstrap interval:** the interval of the depth ratio contains the truth in 38 of 113 (33.6%; interval over isolates 24.3 to 43.6%): 29 of 81
  one-copy units and 9 of 32 units with two or more copies. The interval is narrower than the error of the
  ratio.
- **Graph contexts:** the count equals the truth in 69 of 113 units, is a correct lower bound in 32, is
  below the truth and not reported as a bound in 6, is above the truth in 0, and gives neither count nor bound in 6.
- **Consensus:** status `SUPPORTED` in 89 units, correct in 83; `REVIEW` in 14, correct in 6; `LOWER_BOUND` in
  10, correct in 10 (a continuous dosage counts as correct when it rounds to the truth). The 6 units
  that are `SUPPORTED` and wrong (`copy_number_supported_but_wrong.tsv`): FDAARGOS_1 23S (truth 7, reported 2, depth ratio 4.80); FDAARGOS_149 23S (truth 6, reported 2, depth ratio 5.40); FDAARGOS_149 16S (truth 6, reported 2, depth ratio 4.97); FDAARGOS_199 16S (truth 6, reported 2, depth ratio 5.96); FDAARGOS_144 23S (truth 7, reported 2, depth ratio 7.45); FDAARGOS_209 16S (truth 4, reported 1, depth ratio 1.60).
  In five of them the graph gives two contexts for an rRNA gene that has 6 or 7 copies while the depth
  ratio is between 4.80 and 7.45.

### Block size of the bootstrap (targeted test of protocol section 3, external units)

`estimate_locus_copy_number()` of `locus_recon/depth_ratio.py` was called with `block` of 100, 250, 500 and
1,000 bp and every other argument at its default (`block_size_sensitivity.tsv`). At 500 bp the function returns
the values of the command for all 128 units. Protocol set, 113 units with a truth copy:

| Block (bp) | Median width of the interval | Interval contains the truth |
|---|---|---|
| 100 | 0.130 | 32 of 113 (28.3%) |
| 250 | 0.145 | 37 of 113 (32.7%) |
| 500 | 0.153 | 38 of 113 (33.6%) |
| 1,000 | 0.130 | 30 of 113 (26.6%) |

The block size changes the width little and does not make the interval contain the truth more often.

### Supplementary set (12 units, bait by fallback rule)

12 units with a truth copy: Locus-Recon accepted 8, all exact, and did not accept 4, none exact (the
four *tuf* units with two different copies). KMA accepted the same 8 under both rules and ARIBA 8
(`decision_tables.tsv`, unit set `supplementary`).

### Reproducibility

- A second, independent run of FDAARGOS_262 with the final runner in an empty working directory gave the same
  584 values (sequence checksums, decisions, tiers, flags, ratios, intervals and copy-number calls):
  `rerun_check_FDAARGOS_262.tsv`.
- The stages ran with successive versions of the runner (`PROVENANCE.json`, `stages`, `script_sha256`). The
  changes after the first tool run were resumable downloads, one record file per isolate, the dnadiff stage and
  the scoring and report code; `DEVIATIONS.md` lists them.

## How to rerun

```
# environment of environment.yml (Locus-Recon 1.0.1, SPAdes 4.2.0, BWA-MEM2 2.2.1, samtools 1.22.1, BLAST+ 2.16.0), plus seqtk
python validation/run_external_evaluation.py --workdir W --stage prepare --seqtk /path/to/seqtk
python validation/run_external_evaluation.py --workdir W --stage assemble --assembly-jobs 1 --assembly-threads 8
python validation/run_external_evaluation.py --workdir W --stage reconstruct
python validation/run_external_evaluation.py --workdir W --stage copy-number
python validation/run_external_evaluation.py --workdir W --stage block-size
# environment with ARIBA 2.14.7, KMA 1.7.1 and MUMmer
python validation/run_external_evaluation.py --workdir W --stage comparators
python validation/run_external_evaluation.py --workdir W --stage concordance
# first environment again
python validation/run_external_evaluation.py --workdir W --stage score --outdir validation/external-evaluation-results
```

Each stage skips what is already done. `--isolates` and `--loci` restrict a run. `prepare --verify-subsample`
subsamples a second time and requires identical files. The working directory needs about 55 GB. `host_scripts/`
holds the shell scripts that ran these stages on the host, with its paths and conda environments.

## Choices that the protocol leaves open, fixed before any tool was run on the reads

- Subsampling: `seqtk sample -s 20261118 FILE N` on each mate file, with N = ceil(100 x genome length x pairs /
  bases); the names of the two files were checked to be in step.
- Draft: `contigs.fasta` of SPAdes `--isolate`, as `_draft()` of `run_tier_calibration.py`.
- Truth: `blastn -task blastn -evalue 1e-10` of the bait against the complete genome. An alignment of at least
  80% identity over at least 90% of the bait is one copy, extended base for base to the span of the bait; two
  alignments that overlap by more than half on the same replicon are one copy. With megablast the number of
  copies is the same in all 128 units (`baits.tsv`, `truth_copies_megablast`).
- Intervals for the copy-number commands: from `blast_step1.tsv` of the Locus-Recon run, with the selection of
  `_discover()` of `validation/workflows/run_low_copy_23S.py` (identity of at least 80%, length of at least
  max(200 bp, 10% of the bait), highest bit score first, no overlap on a contig). `--gene-span` and
  `--gene-length` were always passed; `--reference` was not.
- A continuous dosage is scored as correct when it rounds to the true copy number.

## Files

| File | Content |
|---|---|
| `external_units.tsv` | One row per unit and tool (384 rows): truth stratum, returned sequence, relation to the truth, decision under each rule, and for Locus-Recon the disposition, tier and flags |
| `decision_tables.tsv` | Decision by result for each rule, stratum and unit set |
| `external_summary.tsv` | Every proportion with its interval over isolates, and medians (656 rows) |
| `external_headline.tsv` | The rows of the summary that answer protocol 2.4 to 2.6 (261 rows) |
| `counts_per_isolate.tsv`, `counts_per_locus.tsv` | Counts per isolate and per locus for each rule |
| `accepted_not_identical.tsv` | Sequences accepted by Locus-Recon that differ from the truth |
| `mixed_sites.tsv` | Every credible mixed site of every candidate, with the bases of the truth copies at the site |
| `locus_recon_report_rows.tsv` | The full report row of every Locus-Recon run |
| `external_copy_number.tsv` | One row per unit: depth ratio, interval, call, dosage, graph contexts, consensus, and each against the truth |
| `copy_number_error_by_factor.tsv` | Relative error of ratio and dosage by G+C of locus and genome, locus length and number of intervals |
| `copy_number_supported_but_wrong.tsv` | Units whose consensus is `SUPPORTED` and wrong |
| `block_size_sensitivity.tsv` | Ratio and interval at the four block sizes, per unit |
| `baits.tsv`, `baits.fasta` | Source, rule and checksum of every bait; the sequences |
| `truth_copies.tsv`, `truth_sequences.fasta` | Position, strand and checksum of every truth copy; the sequences |
| `locus_recon_candidates.fasta` | The candidate of every Locus-Recon run |
| `comparator_returned_sequences.fasta` | The sequences scored for ARIBA and KMA |
| `input_reads.tsv` | ENA files with MD5 and SHA-256, read and base counts, pairs kept, checksums of the subsampled files |
| `isolates_and_drafts.tsv` | Genome and draft statistics, SPAdes command, and the comparison of draft and genome |
| `rerun_check_FDAARGOS_262.tsv` | Values of the second run beside those of the main run |
| `not_run.tsv` | Anything that did not run (empty) |
| `per_unit/<isolate>/<locus>/` | Candidate FASTA, QC report, `allele_remap.support.tsv.gz`, report row, `blast_step1.tsv`, depth-ratio TSV, copy-number TSV and JSON, block-size JSON, commands |
| `comparator_raw/<isolate>/` | Catalogue, ARIBA report and sequences, KMA `.res` and `.fsa`, logs and commands |
| `concordance_reports/` | dnadiff report of each draft against its complete genome |
| `host_scripts/`, `logs/` | Shell scripts and logs of the run on the host; bait-symbol and tool-version probes |
| `PROVENANCE.json` | Versions, package digest, commands, seeds, parameters, input checksums, stage records |
| `DEVIATIONS.md` | Departures from the protocol |

## Limits

- The truth is not independent of short reads: FDA-ARGOS assemblies combine PacBio and Illumina data.
- All libraries come from one consortium, and one of the sixteen does not match its genome.
- Four other drafts hold sequence that the complete genome lacks: FDAARGOS_780 (97.2% of the draft aligned), FDAARGOS_156 (97.2% of the draft aligned), FDAARGOS_750 (91.8% of the draft aligned), FDAARGOS_501 (84.5% of the draft aligned).
  Its origin was not examined.
- The loci were chosen because every species has them. The units with several different copies are rRNA genes
  and *tuf*; no long or repeated virulence gene is represented, and no amplified resistance gene.
- A copy that spans the origin of a circular replicon, or that an insertion splits, is not found by the truth
  rule.
