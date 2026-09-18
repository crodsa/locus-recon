# Frame, bait validation, placeholder junction and path scoring evidence

This directory validates four reporting behaviours and one ranking behaviour
against material that is already deposited in this repository: the 111 curated
`tcdB` alleles and the LIBA-6656 local assembly graph under
[`../liba6656-gfa-rerun/`](../liba6656-gfa-rerun/README.md). No new sequence
data is introduced. Every derived input is constructed from those two files at
run time and is not stored, except the simulated reads of the positive control,
which are small and deterministic.

## Inputs

- `tcdB_111_ungapped_baits.fasta`: 111 curated alleles, 7,101-7,104 bp.
- `assembly_graph_after_simplification.gfa`: the LIBA-6656 local graph, from
  which 512 terminal paths of 8,000-9,000 bp are enumerated.
- ENA run `ERR467623` for the published scoring case only, verified against the
  MD5 values the ENA file report supplies. The reads are not duplicated here.

`INPUT_CHECKSUMS.tsv` binds the deposited inputs to the files used, and
`PROVENANCE.json` records the commit, parameters and tool versions.

## Command

```bash
python validation/run_frame_and_graph_validation.py \
  --deposit validation/liba6656-gfa-rerun \
  --outdir validation/frame-and-graph-evidence \
  --reads-r1 ERR467623_1.fastq.gz \
  --reads-r2 ERR467623_2.fastq.gz \
  --threads 24
```

Checks 1-5 need no reads and run in seconds; check 6 needs the ENA pair and an
aligner. Omitting the read arguments runs everything else.

## Results

### 1. The frame metrics do not depend on bait orientation

One allele (`B1.1`) is held out as the candidate and the remaining 110 build
the profile. The pair is assessed in all four orientation combinations
(`frame_orientation_invariance.tsv`).

| quantity | result |
|---|---|
| Internal stops, six-frame, all four combinations | **0** |
| Internal stops, forward frames only | **0** when the candidate is on the deposited strand, **103** when it is reverse complemented |
| Frame reported | `+1` and `-1`, mirroring the candidate's orientation |
| Reverse-complemented bait set | reported as reading coding only on the reverse strand |

The 103 stops are the quantity a three-frame scan produces for a candidate that
is a clean open reading frame on the other strand. They are a property of the
scan, not of the sequence, which is why the count is now minimised over six
frames and the frame used is reported.

### 2. Bait database validation names the record and the fix

Three copies of the deposited bait file are derived and submitted
(`bait_validation_messages.tsv`).

| derived input | outcome |
|---|---|
| A 64-character identifier | rejected, naming the record, its length and the 50-character `makeblastdb -parse_seqids` limit |
| Twenty alignment gaps restored in one record | rejected, naming the symbol and that the file looks like a gapped alignment |
| A 400 bp unrelated fragment added to six alleles | accepted with a warning naming the record as a probable paralogue |

The deposited file itself is accepted unchanged, which is the false-positive
control for these three checks.

### 3. Length flags carry the size of the profile behind them

The same truncated candidate is assessed against a three-allele subset and
against all 110 alleles (`length_profile_reporting.tsv`). Both report
`LENGTH_ANOMALOUS (delta=3551 bp, outside 3xIQR)`; only the three-allele
profile adds `profile n=3` with
`CATALOGUE_LENGTH_PROFILE_UNDERPOWERED`, and the 110-allele profile reports
`profile n=110` without it. The flag text is the only difference, which is the
point: the tier is unchanged and the reader can see how well determined the
comparison was.

### 4. Graph evidence for a scaffold placeholder, three ways

512 terminal paths are enumerated from the deposited graph; the longest is
8,538 bp. Three candidates are constructed from it and assessed
(`placeholder_junction_cases.tsv`).

| case | expected | observed | graph paths carrying both flanks | distance between flanks |
|---|---|---|---|---|
| 100 bp placeholder inserted where the graph is contiguous | `GRAPH_SUPPORTED` | `GRAPH_SUPPORTED` | 256 | 0 bp |
| 100 bp placeholder replacing 150 bp the graph still spells | `AMBIGUOUS` | `AMBIGUOUS` | 128 | 150 bp |
| Second flank shuffled so it is absent from the graph | `NOT_SUPPORTED` | `NOT_SUPPORTED` | 0 | none |

Three of three cases agree. The first case is the one a local assembler
produces when it scaffolds across a gap the graph itself closes; the second is
the case where the graph offers a different spelling and the junction is
therefore not evidence; the third is the negative control.

### 5. Path ranking recovers a known source path

The published case in check 6 cannot separate its paths, so the ranking is also
exercised against constructed truth, in the same spirit as the completeness
benchmark. Reads are simulated from one deposited path (40x, 2 x 100 bp, 0.2%
substitutions, seed 7, 1,707 pairs) and scored against that path plus the nine
most divergent paths in the deposit
(`scoring_positive_control.tsv`).

| candidate | mapped reads | uniquely anchored depth | uniquely anchored breadth |
|---|---|---|---|
| source (rank 1) | 636 | 3.42x | **10.21%** |
| decoy_1 to decoy_9 (ranks 2-10) | 255-367 | 0.00x | **0.00%** |

The source path ranks first and is the only candidate with any uniquely
anchored coverage. Its breadth is 10.21% rather than a high number because the
decoys are alternative paths through the same graph and share the remaining
90% of its sequence: the metric measures the sequence that distinguishes a
path, not the sequence the candidates have in common. Every decoy still
attracts 255-367 reads, which is why a mapped-read count alone cannot rank
these candidates.

### 6. On the published case the reads cannot separate the paths

All 512 enumerated paths were scored against `ERR467623`
(`graph_path_scores.tsv`).

| quantity | result |
|---|---|
| Paths scored | 512 |
| Mapped reads per path | 10-43 |
| Uniquely anchored breadth | **0.00% for every path** |
| Verdict reported by the command | the paths cannot be told apart at this read length, and their order is arbitrary |

Every read fits several paths equally well and is therefore excluded by the
mapping-quality filter. This is the expected result for a locus whose
alternatives differ only at graph bubbles shorter than the library insert, and
it is the measured form of the boundary the original deposit states in prose:
each exported path is a hypothesis. The ranking adds an ordering where one
exists and says so when none does.

## Boundaries

- Checks 1-5 are deterministic and reproduce byte for byte from the deposited
  files; check 6 depends on the ENA pair remaining available.
- The positive control uses simulated reads: it establishes that the ranking
  recovers a known source among divergent alternatives, not a sensitivity for
  real libraries.
- Check 4 constructs the placeholders rather than observing one, so it
  validates the junction logic against the real graph, not the frequency with
  which a local assembler inserts such a placeholder.
- None of these checks change the confidence tiers. The placeholder and length
  additions are descriptive by construction, and the two frame corrections
  remove flags that described the scan rather than the sequence.
