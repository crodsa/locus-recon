# Contemporary LIBA-6656 GFA rerun

This directory contains a contemporary reconstruction of the local SPAdes graph for the LIBA-6656 `tcdB` analysis. It was generated on 24 August 2026 after the original paired reads were recovered. It is not represented as the exact historical GFA used during the initial exploratory analysis.

## Inputs

- ENA run `ERR467623`: 1,383,413 synchronized 2 x 100 bp read pairs.
- Draft assembly `LIBA6656_ST154.fasta`: 77 contigs and 4,268,222 bp.
- Bait source `tcdb_diffbase.curated_alignment.fasta.gz`: 111 curated sequences. The deposited bait FASTA was produced only by deleting alignment-gap characters (`-`); sequence identifiers and all nucleotide characters were retained.

The MD5 and SHA-256 values in `INPUT_CHECKSUMS.tsv` bind these inputs to the files used in the rerun. Raw reads are not duplicated here because they are available through ENA accession `ERR467623`.

## Software and command

The local assembly used Locus-Recon 1.0 (commit `9c02bd1`), BLAST+ 2.16.0+, BWA-MEM2 2.2.1, samtools 1.22.1, and SPAdes 4.2.0. Discovery used 80% minimum identity, a 50 bp minimum HSP, and 300 bp interval padding. Eight threads and 16 GB SPAdes memory were requested. Cleanup was disabled.

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

## Files

- `assembly_graph_after_simplification.gfa`: contemporary SPAdes graph.
- `tcdB_111_ungapped_baits.fasta`: exact ungapped bait input.
- `tcdB_terminal_paths_8000_9000.fasta`: all 512 terminal paths.
- `tcdB_terminal_paths_8000_9000.tsv`: path lengths, graph depths, and oriented node lists.
- `selected_terminal_paths_210_315.fasta`: the two paths described above.
- `reconstructed_candidates_vs_terminal_paths.tsv`: BLAST alignments of both deposited candidates against every path.
- `GFA_RERUN_PROVENANCE.json`: portable machine-readable parameters and outcome summary.
- `INPUT_CHECKSUMS.tsv`: checksums for source inputs and deposited rerun outputs.
