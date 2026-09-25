# Pending numbers before the 1.0.0 tag

The documentation carries slot tokens, `[[CALIB-nn]]` or `[[CALIB-nn note]]`,
wherever a number from tier-calibration stages B to D belongs. Those stages run
on public *H. pylori* reads and RefSeq annotations and must be re-run with the
1.0.0 code: the per-base pileup now counts reads that are not properly paired
(`samtools mpileup -A`), the tier rule is monotone and scores a contig-end
truncation once, and 23S is reconstructed with `--noncoding-locus`. No slot may
be filled from the previous run's tables.

## 1. Run the stages

From the repository root, in the `locus-recon` environment (add `matplotlib`
if the figure `tier_calibration.{pdf,svg}` is wanted; without it the run
completes and draws no figure):

```bash
python validation/run_tier_calibration.py \
  --outdir validation/tier-calibration \
  --manifest validation/tier-calibration/closed_genome_manifest.tsv \
  --workdir /scratch/tier-calibration --threads 24 \
  --divergent-loci 3 --depth-series 15,8 --depth-isolate
```

The run needs NCBI (RefSeq assemblies and annotations) and the ENA or SRA read
archives. It rewrites the stage tables, `TIER_CALIBRATION.json`,
`tier_calibration_all_cases.tsv`, `tier_calibration_tables.md`,
`PROVENANCE.json` and `INPUT_CHECKSUMS.tsv`.

## 2. Close the slots

```bash
python release/slot_gate.py --values   # prints the value of every slot
python release/slot_gate.py --list     # where each slot occurs
```

Write each value in place of every occurrence of its token, and reread the
sentence around it: a slot that comes back as zero, or larger than the sentence
assumes, means the sentence has to change, not only the number.

| Slot | Quantity | What closes it | Kind |
|---|---|---|---|
| `[[CALIB-01]]` | Exact reconstructions among the ten stage-B loci (23S and *gyrB*) | `--values`, from `tier_calibration_closed_genomes.tsv` | compute |
| `[[CALIB-02]]` | *gyrB* reconstructions exact and at `HIGH`, of five | `--values`, stage B table | compute |
| `[[CALIB-03]]` | Identity range of those *gyrB* reconstructions to the catalogue allele | `--values`, `catalogue identity` in the stage B evidence | compute |
| `[[CALIB-04]]` | Stage-B reconstructions that were not exact: case, tier and disposition, relation to the annotation, completeness flags | `--values`; phrase it as prose in README.md and the calibration README | compute |
| `[[CALIB-05]]` | Supported cases at `HIGH`, of supported cases, all four stages | `--values`, `combined` in `TIER_CALIBRATION.json` | compute |
| `[[CALIB-06]]` | Precision of `HIGH` with its Clopper-Pearson 95% interval | `--values`, `combined.overall` | compute |
| `[[CALIB-07]]` | False accepts, all four stages | `--values`, `combined.overall` | compute |
| `[[CALIB-08]]` | Withheld cases exact over the span they reported, of withheld cases | `--values`, `combined.overall` | compute |
| `[[CALIB-09]]` | Reported sequences whose bases disagree with truth, with depth, tier and flags | `--values`, `relation_to_annotation` of stages B to D | compute |
| `[[CALIB-10]]` | Depth series: tiers and exact sequences at full depth, ~15x and ~8x | `--values` prints the table | compute |
| `[[CALIB-11]]` | Divergence arm: per-locus identity to the bait, bait length, cases at `HIGH`, relation to truth | `--values` prints the table | compute |
| `[[CALIB-12]]` | Two-copy 23S loci: tiers, and consensus exact to the annotation | `--values`, stage B 23S rows | compute |
| `[[CALIB-13]]` | Per-class table for all four stages | `--values`, or paste "All four stages" from `tier_calibration_tables.md` | compute |
| `[[CALIB-14]]` | Tier agreement of the `--isolate` control with the primary depth series | `--values`, `stage_d_isolate_control` | compute |

All fourteen close by computation; none is an author decision or depends on
data outside the run.

## 3. Propagation that no slot carries

These items are invisible to the gate and are the ones that get forgotten.

- **Wording that depends on the values.** README.md: the closed-genome bullet
  near the top ("which is the case the separate catalogue axis exists for")
  holds only if `[[CALIB-02]]` is above zero; the "What a tier is worth"
  statements and the tier-calibration README statements 1 to 4 frame
  questions whose answers are the slots. The Stage C text and the class
  assignment (`divergent` below 95% identity) are rules; check that the
  loci the run selects are described correctly.
- **Low-copy 23S deposit.** `validation/low-copy-23S/` was produced before the
  released consensus labels. Under 1.0.0 the five 23S rows, whose graph count
  of two lies outside their depth interval, carry
  `copy_number_method = GRAPH_COUNT_OVER_DISCORDANT_DEPTH` rather than
  `GRAPH_DEPTH_CONSENSUS`, and the table carries two descriptive columns
  (`profile_call`, `profile_calibrated`) that 1.0.0 no longer writes. Every
  criterion is unaffected. It uses the same five genomes as stage B:
  `python validation/workflows/run_low_copy_23S.py --work-dir work/low-copy-23S --threads 4 --memory-gb 8`,
  then copy the compact outputs into `validation/low-copy-23S/` as before. If
  that re-run is not possible, add an editorial note to
  `validation/low-copy-23S/PRESPECIFICATION.md` stating both differences.
- **LIBA-6656 reconstruction statement.** README.md ("the standard workflow
  rejected a 1,955 bp mixed and frame-anomalous *tcdB* sequence") comes from
  the reconstruction run behind `validation/liba6656-gfa-rerun/`. The sequence,
  the frame anomaly and the rejection do not depend on the pileup change, but
  the mixture call does; re-running that reconstruction with 1.0.0 on
  `ERR467623` and the archived draft `LIBA6656_ST154.fasta` confirms it.
- **Dates.** `CITATION.cff` (`date-released`) and the `CHANGELOG.md` heading
  carry 2026-09-25; set both to the day the final tag is made.
- **The tag.** `v1.0.0` points at a commit that still carries slots, and the
  `release-gate` workflow fails on it by design. After closing the slots:
  `git tag -f v1.0.0 && git push -f origin v1.0.0`.
- **This directory.** Delete `release/SLOTS.md` once every slot is closed;
  keep `release/slot_gate.py` and `.github/workflows/release-gate.yml` for
  later releases, or delete all three together.

## 4. Gate

`python release/slot_gate.py` exits non-zero while any slot remains in the
README, CHANGELOG, CITATION, `docs/`, `examples/` or `validation/` Markdown.
The `release-gate` workflow runs it on every `v*` tag.
