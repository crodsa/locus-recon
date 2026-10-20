# Contemporary LIBA-6656 GFA rerun

This directory contains a contemporary reconstruction of the local SPAdes graph for the LIBA-6656 `tcdB` analysis. It was generated on 24 August 2026 after the original paired reads were recovered. It is not represented as the exact historical GFA used during the initial exploratory analysis.

## Inputs

- ENA run `ERR467623`: 1,383,413 synchronized 2 x 100 bp read pairs.
- Draft assembly `LIBA6656_ST154.fasta`: 77 contigs and 4,268,222 bp.
- Bait source `tcdb_diffbase.curated_alignment.fasta.gz`: 111 curated sequences. The deposited bait FASTA was produced only by deleting alignment-gap characters (`-`); sequence identifiers and all nucleotide characters were retained.

The MD5 and SHA-256 values in `INPUT_CHECKSUMS.tsv` bind these inputs to the files used in the rerun. Raw reads are not duplicated here because they are available through ENA accession `ERR467623`.

## Software and command

The local assembly used Locus-Recon 1.0.0, BLAST+ 2.16.0+, BWA-MEM2 2.2.1, samtools 1.22.1, and SPAdes 4.2.0. Discovery used 80% minimum identity, a 50 bp minimum HSP, and 300 bp interval padding. Eight threads and 16 GB SPAdes memory were requested. Cleanup was disabled.

```bash
python -m locus_recon.cli \
  --samplesheet samples.tsv \
  --main-output-dir output \
  --bait tcdB_111_ungapped_baits.fasta \
  --locus tcdB \
  --threads 8 \
  --parallel-samples 1 \
  --memory-per-sample 16 \
  --min-discovery-identity 80 \
  --min-hsp-length 50 \
  --region-padding 300 \
  --no-progress
```

Terminal graph paths were then exported with:

```bash
python -m locus_recon.graph_cli \
  --gfa assembly_graph_after_simplification.gfa \
  --output tcdB_terminal_paths_8000_9000.fasta \
  --summary tcdB_terminal_paths_8000_9000.tsv \
  --min-length 8000 \
  --max-length 9000 \
  --prefix LIBA6656_tcdB_terminal_path
```

## Result and interpretation

The GFA contains 26 segments and 32 links. The exporter retained all 512 terminal paths, with lengths from 8,456 to 8,538 bp. The count `512 = 2^9` is consistent with the nine variable graph regions described in the article. These 512 combinatorial paths are unresolved proposals, not 512 biological alleles.

Alignment of the deposited 7,104 bp reconstructed candidates against all paths identified:

- `LIBA6656_tcdB_terminal_path_210`: exact full-length match to the extrachromosomal-context candidate, in reverse orientation.
- `LIBA6656_tcdB_terminal_path_315`: best full-length graph match to the chromosome-associated candidate (99.226% identity; 55 differences). Those 55 differences correspond to the subsequent polishing stage and are not silently written into the GFA path.

The two selected terminal paths are provided for navigation only. Candidate acceptance still rests on the deposited competitive-remapping, per-base, context, and polishing evidence. The contemporary graph therefore restores an auditable graph/path enumeration while preserving the distinction from the unavailable historical graph file.

## Verification with the release code

On 28 September 2026 the two commands above were rerun on the release code
(source digest `3863abc32fb5446a72919c790cb2558217871410f441db1ecf17dde3de55cdcf`)
from inputs whose checksums equal those in `INPUT_CHECKSUMS.tsv`, once with four
threads and once with eight, with a 12 GB SPAdes memory cap, under Python
3.11.13 and the tool versions listed above. `verify_release_rerun.py` compares
the reruns with this deposit and writes `RELEASE_VERIFICATION.json`.

- Both reruns recruited the same 6,195 read names and wrote byte-identical
  graphs and reported sequences.
- Their graph is the deposited graph up to numbering: the same 26 segment
  sequences and the same 32 links once segments are matched by sequence.
  Segment identifiers differ, and two segments differ in k-mer coverage, by 2
  of 27,945 and by 8 of 69,108 k-mers.
- The 512 terminal paths exported from it are identical to the deposited paths
  in sequence and numbering, so paths 210 and 315 and the candidate alignments
  above are unchanged. In the path summary only the node identifiers differ,
  and the mean graph depth by at most 0.001.
- The standard workflow reported one 1,955 bp `tcdB` sequence at `SUSPECT`,
  disposition `HOLD`. It carried `ALLELE_MIXTURE` (23 bidirectional sites,
  median alternative fraction 0.404) and `ALLELE_SPAN_CLIPPED_AT_CONTIG_END`
  (1,345 bp at the contig start and 3,804 bp at its end, 5,149 bp in all);
  `LOCUS_SPLIT_ACROSS_CONTIGS` located 1,345 of the missing bases on a second
  contig, and `FRAME_LENGTH_SHIFT_TRUNCATED` attributes its length, which is
  not a multiple of three, to that truncation. Mixed sites near an allele
  fraction of 0.4 are consistent with two divergent copies collapsed into one
  contig, the situation the two deposited candidates resolve; the workflow
  withholds the consensus rather than reporting it as an allele.

## Files

- `assembly_graph_after_simplification.gfa`: contemporary SPAdes graph.
- `tcdB_111_ungapped_baits.fasta`: exact ungapped bait input.
- `tcdB_terminal_paths_8000_9000.fasta`: all 512 terminal paths.
- `tcdB_terminal_paths_8000_9000.tsv`: path lengths, graph depths, and oriented node lists.
- `selected_terminal_paths_210_315.fasta`: the two paths described above.
- `reconstructed_candidates_vs_terminal_paths.tsv`: BLAST alignments of both deposited candidates against every path.
- `GFA_RERUN_PROVENANCE.json`: portable machine-readable parameters and outcome summary.
- `INPUT_CHECKSUMS.tsv`: checksums for source inputs and deposited rerun outputs.
- `standard_workflow_report_tcdB.tsv`, `standard_workflow_qc_report_tcdB.txt`,
  `standard_workflow_reconstructed_tcdB.fasta`: the report row, QC report and
  withheld sequence of the release rerun; `allele_file` is shown relative to
  the results directory.
- `verify_release_rerun.py`, `RELEASE_VERIFICATION.json`: the comparison of the
  release reruns with this deposit, and its record.
