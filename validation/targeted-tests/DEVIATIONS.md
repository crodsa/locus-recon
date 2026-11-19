# Departures from the protocol: targeted tests with simulated reads

Protocol: `validation/external-evaluation/PROTOCOL.md`, section 3 (written on 18 November 2026).
Runner: `validation/run_targeted_tests.py`. Frozen software: Locus-Recon 1.0.1, package digest
`2cc816cca467a7988c1f2a28cba11ed1cbb30b263e331f75ca48efa19ebefce0` (checked by the runner before
any run, and recorded in `PROVENANCE.json`). No file of the package, no threshold and no default
was changed.

Each entry states what was done, what the protocol says, and why. Entries D1 to D7 were fixed
before any sample was simulated. Entries D8 to D11 were decided after three pilots and before the
full run. Entry D12 is the only change made after the full run had started.

## Design choices that the protocol leaves open or that depart from its wording

**D1. Read sampling.** The protocol says that reads carry substitution errors "as in
`validation/generate_mock_dataset.py`". The error model of the generator is kept unchanged
(150-bp pairs, substitution rate 0.001, base quality `I`). Fragment sampling is extended:
fragment starts are uniformly random (as in `validation/workflows/run_constructed_copy_series.py`),
the orientation of each pair is random, and fragment length is normal with mean 350 bp and
standard deviation 30 bp, limited to 200 to 500 bp. Reason: the generator places fragments at
even spacing with alternating orientation, so replicates would differ only by a jitter of 8 bp
and would not sample depth noise or strand noise, which are the object of tests 1 and 2.

**D2. Draft assemblies.** The repository's mock benchmark gives Locus-Recon drafts cut from the
true genome. Here every simulated sample is assembled with `spades.py --isolate` (SPAdes 4.2.0,
one thread), as `validation/workflows/run_low_copy_23S.py` does, and that draft and its graph
are the inputs. The exceptions are the drafts of test 5, which the protocol itself defines as
"cut" (see D10 and D11).

**D3. Locus lengths.** The alleles are coding sequences built as the generator builds them
(`ATG` plus sense codons), so their length is a multiple of three. The 500-bp and 7,000-bp loci
of test 2 are 501 bp and 6,999 bp. Tests 1 and 3 use a 1,500-bp locus and tests 4 and 5 a
3,000-bp locus; the protocol states no length for them. The background is 100,000 bp of random
sequence.

**D4. Bait.** The bait is one sequence. In test 1 it is the major allele; in tests 2, 4 and 5
it is the allele of the sample (copy 1 where copies differ); in test 3 it is the reference
allele from which the sample allele was derived.

**D5. Test 3, scenarios.** "Substitutions within 5, 20 or 50 bp of an end" is implemented as one
substitution within that distance of the 5' end and one within that distance of the 3' end of
the same allele. A ninth scenario, an allele identical to the bait, was added as a control
(15 runs). Indels and the stop codon are placed at least 200 bp from either end. The
compensating indels are a 1-bp deletion and a 1-bp insertion 30 bp downstream. Whether they
leave a stop codon in the ten altered codons depends on the replicate; the truth of each run is
in the table (`truth_frame_intact`, `truth_internal_stops_frame_plus1`).

**D6. Discovery intervals for the copy-number commands.** `locus-recon-depth-ratio` and
`locus-recon-copy-number` need `--locus` and `--gene-span`. They are taken from a BLASTN of the
bait against the draft with the rule of `validation/workflows/run_low_copy_23S.py`
(identity of at least 80%, hits of at least max(200 bp, 10% of the bait), best hit kept where
hits overlap on the subject). The alignment is `<sample>/full_map.sorted.bam` written by
`locus-recon`.

**D7. Consensus with and without the graph.** `locus-recon-copy-number` was run twice on every
sample of tests 4 and 5: with the SPAdes graph of the draft (`--gfa`, `--bait`) and without it.
The second call is the configuration of the *tcdB* case of LIBA-6656. Both are reported
(`cn_*` and `cn_nograph_*` columns).

## Changes decided after the pilots, before the full run

Three pilots were run to test the runner. They used replicate number 91, whose seeds are not
used by any run of the evaluation (replicates 1 to 20). Their tables are in `pilot/`. They are
not part of the results. The pilots led to these changes.

**D8. Second background for test 4 (decoy repeat).** In the design written before the pilot the
locus was the only repeat of a linear genome. In the first pilot SPAdes then returned one
scaffold that holds both copies of a two-copy locus, also for a tandem duplicate, because the
graph of such a genome has a single traversal. This does not happen in a bacterial genome, which
has other repeats. A second background was therefore added: an unrelated 2,000-bp repeat in two
copies, placed at the midpoint of the last two background segments so that it interleaves with
the locus copies. Test 4 was run in full on both backgrounds (`background` = `plain`, the design
as first written, run identifiers `t4_`; and `decoy`, run identifiers `t4d_`): 170 runs
in place of 85. Both are reported and neither is given precedence.

**D9. Test 5 uses the decoy background.** The cut that the protocol asks for needs a contig that
holds the whole locus once. With the plain background the two-copy draft holds two complete
copies (pilot 1: 2 of 2 two-copy samples could not be cut), so both copy numbers of test 5 use
the decoy background, in which the two copies collapse into one contig.

**D10. Test 5, how the draft is cut.** The contig that holds the whole locus is replaced by two
pieces: piece A ends at query position *b*, piece B starts at query position *a*, and
*b* - *a* + 1 = round(overlap x 3,000), centred on the query. At 0% the two pieces are adjacent.
For two copies the contig is the collapsed locus contig. The overlapped bases are therefore
present twice in the draft with identical sequence.

**D11. Test 5, added arm with two divergent copies.** Not in the protocol. In the arms of D10
the reads of the overlapped part map equally well to both pieces. In the *tcdB* case that
motivated the test the overlapping intervals come from two copies that differ, and reads map
uniquely. An arm was added to cover that geometry with known truth: two copies that differ at
1% of sites; the draft keeps copy 1 over query 1 to *b* and copy 2 over query *a* to 3,000 and
lacks the rest of both (35 runs, `arm` = `c2_divergent_partial`, run identifiers `t5_c2div_`).
SPAdes returned both divergent copies complete in only 2 of 4 pilot samples, so the draft of
this arm is built from the true chromosome, as the drafts of the repository's mock benchmark
are; SPAdes was still run on the reads and its graph is the one given to
`locus-recon-copy-number`. The two arms that the protocol defines are unchanged and are reported
in full.

## Change made after the full run had started

**D12. Four runs of the added arm were repeated after a correction of the runner.** In the full
run, 4 of the 35 runs of the added arm (D11) ended as `NOT_CUT_NEEDS_TWO_COMPLETE_COPIES`
(`t5_c2div_ov10_r05`, `t5_c2div_ov65_r01`, `t5_c2div_ov65_r04`, `t5_c2div_ov90_r02`). The cause
was in the runner and not in Locus-Recon: a substitution within three bases of an end of copy 2
shortened the BLASTN alignment that the cut function uses to find copy 2, and the function
required a full-length alignment. The function now extends such an alignment to the whole query
(the copies differ by substitutions only). These four runs were repeated with the same seeds;
Locus-Recon had not been run on them before. Their first records are kept on the host
(`work/superseded_not_cut/`). The other 786 runs were not repeated: they were produced by the
first version of the runner (SHA-256 `371175a4bf0c9ef603923aa6de4084c15db7e170c6908883ad9ce4f32a1e890f`),
and `logs/runner_changes_after_full_run.diff` lists every difference between that version and
the delivered one. The differences are the cut function of the added arm and the summary tables
(marginal totals); the simulation, the commands and the scoring of the other tests are the same.

## Other notes

**D13. Progress display.** `locus-recon` was run with its default options plus `--threads 1` and
`--memory-per-sample 8`. `--no-progress` was not used.

**D14. Parts of section 3 that belong to other tracks.** The strand-bias rate in "every one-copy
unit of the external evaluation" and the block-size series on the *aphA1* runs, the 23S rRNA
gene of *H. pylori* and the external units are not in this folder. This folder adds block sizes
of 100, 250, 500 and 1,000 bp for the simulated samples of tests 4 and 5 only.

**D15. Where the work was run.** All runs were made on one Linux host with eight runs in
parallel, one thread each (`logs/full_run_time.txt`). Local execution in the analysis sandbox
was not available during this work, so no part of the runner was executed outside that host.

**D16. Limits of the simulation that affect how the results can be read.** Errors are
substitutions at a uniform rate, equal on both strands, with constant base quality. Test 2
therefore measures strand-bias alerts that arise from sampling alone; it does not contain
strand-specific or context-specific errors. Genomes are linear random sequences of 100 kb with
50% G+C. In test 4 and test 5, SPAdes resolved some repeats that it would not resolve in a
bacterial genome; D8 and D9 describe how this was handled.
