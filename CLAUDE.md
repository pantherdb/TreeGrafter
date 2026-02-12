# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TreeGrafter is a phylogenetic annotation tool that assigns uncharacterized protein sequences to positions in annotated PANTHER phylogenetic trees. It combines HMMER-based sequence matching with RAxML phylogenetic tree grafting to predict GO terms, subfamily assignments, and protein class annotations.

**Language:** Primarily Perl (`treeGrafter.pl`), with Python utility scripts. No build system — scripts run directly.

## Running

```bash
# Basic run with auto-detection of HMMER algorithm
perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/sample.out -d ./Test/PANTHER_mini -auto

# With pre-computed HMMER results (skips HMMER stage)
perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/sample.out -d ./Test/PANTHER_mini -algo hmmscan -hmmer ./Test/sample.fasta.hmmscan.out

# Docker
docker run --rm -v /path/to/Test:/tmp ningzhithm/treegrafter:1.01 -f /tmp/sample.fasta -o /tmp/sample.out -d /tmp/PANTHER_mini -auto
```

Key flags: `-f` input FASTA, `-o` output, `-d` PANTHER data dir, `-algo hmmscan|hmmsearch`, `-auto` auto-select algorithm, `-hmmer` precomputed results, `-k` keep temp files, `-t` tmpdir, `-cpus` HMMER threads.

## Testing

No test framework. Verify by comparing output against reference:
```bash
perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/test.out -d ./Test/PANTHER_mini -auto
diff ./Test/sample.output ./Test/test.out
```

Test data lives in `Test/` — `sample.fasta` (input), `sample.output` (expected), `sample.fasta.hmmscan.out` (precomputed HMMER).

## External Dependencies

- **RAxML 8.2.4** — phylogenetic inference (`raxmlHPC-SSE3`)
- **HMMER 3.1b2** — HMM sequence searching (`hmmscan`/`hmmsearch`)
- **Perl modules:** Bio::TreeIO, Try::Tiny, JSON::Parse, IO::String
- **Python packages** (for utility scripts only): biopython, numpy

## Architecture — Pipeline Stages

All pipeline logic is in `treeGrafter.pl` (~870 lines), executed sequentially:

1. **HMMER Search** (`runhmmer`) — Runs hmmscan or hmmsearch against PANTHER HMM database. `autodetect` chooses algorithm by comparing fasta size vs HMM database size (threshold: fasta > hmm×40 → hmmsearch).

2. **HMMER Parsing** (`parsehmmer`) — Extracts best-scoring PANTHER family match per query. Different parsing paths for hmmscan vs hmmsearch output formats. Handles multi-domain overlaps (highest-scoring domain wins).

3. **MSA Construction** (`graftMatches` → `_querymsf`) — Aligns each query sequence into its matched PANTHER family's multiple sequence alignment using HMM match states.

4. **RAxML Grafting** (`_runRAxMLAndAnnotate` → `_generateFasta`) — Runs RAxML (`-f y` placement mode) with WAG+GAMMA model to place query into the family's reference tree. Parses jplace (JSON) output.

5. **Annotation Lookup** (`mapto` → `commonancestor` → `printResults`) — Maps graft placement to tree node, finds common ancestor if multiple placements, looks up PANTHER annotations (SF/GO/PC).

## Other Scripts

- `precompute.pl` — Preprocesses PANTHER library data (annotation.dat, node.dat, etc.) into lookup structures
- `cgi-bin/treeGrafter.cgi` — Web service wrapper; reads paths from `config` file (see `config.example`), returns XML
- `sto_to_fasta.py` — Converts Stockholm alignments to FASTA; supports SLURM distribution
- `gettreesandmsf_hpc.pl` — HPC batch processing for PANTHER books
- `fix_bifurcate_files.py` — Strips outer parentheses from Newick tree files

## PANTHER Data Layout

The `-d` data directory expects:
- `famhmm/binHmm` — HMM profile database (with .h3f/.h3i/.h3m/.h3p index files)
- `Tree_MSF/` — Family MSAs as `.AN.fasta` and trees as `.newick`/`.bifurcate.newick`
- `PAINT_Annotations/PAINT_Annotatations_TOTAL.txt` — Annotation mappings (note: typo "Annotatations" is intentional in filename)

## Output Format

Tab-separated: gene_id, PANTHER_family, annotations (semicolon-separated SF/NAME/GOA/PC), graft_point (node ID).

## Further Documentation

See `docs/` for detailed documentation:
- `docs/ARCHITECTURE.md` — In-depth architecture and design details
- `docs/DATA.md` — Data formats, PANTHER library structure, and file specifications
- `docs/TESTING.md` — Testing strategies and verification procedures
