# TreeGrafter Architecture

## Overview

TreeGrafter is a phylogenetic annotation pipeline that assigns functional annotations (GO terms, subfamily, protein class) to uncharacterized protein sequences by placing them into annotated PANTHER phylogenetic trees. The pipeline combines HMMER-based sequence matching with RAxML phylogenetic tree grafting.

The entire core pipeline lives in a single Perl script: `treeGrafter.pl` (~870 lines). There is no module system, no class hierarchy, and no build system. Scripts run directly.

## Component Inventory

| File | Language | Role |
|---|---|---|
| `treeGrafter.pl` | Perl | Core pipeline (HMMER, MSA, RAxML, annotation) |
| `precompute.pl` | Perl | Preprocesses PANTHER library into annotation lookup files |
| `cgi-bin/treeGrafter.cgi` | Perl | Web service wrapper; calls `treeGrafter.pl` |
| `sto_to_fasta.py` | Python | Converts Stockholm alignments to FASTA for PANTHER families |
| `gettreesandmsf_hpc.pl` | Perl | HPC batch extraction of trees and MSAs from PANTHER library |
| `fix_bifurcate_files.py` | Python | Strips outer parentheses from `.bifurcate.newick` files |
| `treeGrafter101` | Binary | Pre-compiled macOS executable (x86_64 + i386 universal binary) |
| `Dockerfile` | Docker | Container build for reproducible execution |
| `config.example` | Config | Template for CGI paths (RAXML_PATH, HMMER_PATH, PTHR_DATA_DIR) |

## Pipeline Stages (treeGrafter.pl)

The main script executes five stages sequentially:

```
[Input FASTA] --> (1) HMMER Search --> (2) Parse HMMER --> (3) MSA Construction
                                                               |
                                                               v
                              [Output TSV] <-- (5) Annotation Lookup <-- (4) RAxML Grafting
```

### Stage 1: HMMER Search (`runhmmer`, lines 210-226)

Runs either `hmmscan` or `hmmsearch` against the PANTHER HMM database.

- **Algorithm selection** is done by `autodetect()` (lines 808-843): reads the total model length from the HMM file (`LENG` lines), multiplies by 40 (approximation of HMM-to-sequence memory ratio), and compares against fasta file size. If fasta size exceeds `hmm_size * 40`, it uses `hmmsearch`; otherwise `hmmscan`.
- Output written to `<fastafile>.<algo>.out`.
- Skipped entirely if the output file already exists with nonzero size or if `-hmmer` flag provides precomputed results.

**Discussion:**
- The autodetect heuristic (40x multiplier) is undocumented in bioinformatics literature. It appears to be an empirical rule of thumb. The threshold may not generalize well across datasets of different composition.
- The `runhmmer` function silently skips execution if output file exists. This is a caching mechanism but has no invalidation -- stale results from a previous run with different parameters will be reused without warning.
- The function redirects stdout to `/dev/null` but not stderr. Errors from HMMER itself may be silently swallowed depending on context.

### Stage 2: HMMER Parsing (`parsehmmer`, lines 232-383)

Extracts best-scoring PANTHER family match per query sequence. Two completely separate parsing paths for hmmscan vs hmmsearch output formats.

**hmmscan path (lines 239-289):**
- Iterates line-by-line through HMMER output.
- Uses regex to identify query sequences (`^Query:`) and PANTHER hits (`^>> PTHR`).
- Only considers the first PANTHER hit per query (`$matchn==1`).
- Extracts HMM coordinates, alignment strings, and domain scores by positional splitting.
- Stores multi-domain data as parallel arrays (hmmstart, hmmend, score, hmmalign, matchalign).

**hmmsearch path (lines 291-382):**
- Two-pass approach: first pass collects all sequence-to-family scores, picks the best family per sequence; second pass extracts alignment details for the winning family.
- Uses hash `%tmpstore{sequence}{family} = score` to find the top-scoring family per sequence.

**Discussion:**
- The two parsing paths share significant duplicated logic (lines 354-380 nearly identical to 261-288) but are maintained separately -- a maintenance burden and source of subtle divergence.
- **hmmsearch path does not extract `score` per domain** (lines 358-360 lack `push(@{...->{score}}, ...)`), but `_querymsf` (Stage 3) sorts by score. This means hmmsearch results will fail or produce undefined behavior when multi-domain overlap resolution is attempted.
- The parser relies heavily on positional column splitting (`split(/ +/)`) of HMMER's human-readable output rather than using HMMER's tabular output formats (`--tblout`, `--domtblout`). This is fragile -- HMMER output formatting can vary between versions.
- Line counter `$hmmline` tracks relative position within each domain block but is not reset between domains within the same query-hit pair. The `$hmmline >= 4` check (line 262) assumes a fixed header size that could break with HMMER version changes.
- The `$matchn==1` guard means only the top-scoring PANTHER family match is considered per query, ignoring potential gene fusion events where a sequence legitimately matches multiple families. A TODO comment at line 259 acknowledges this.

### Stage 3: MSA Construction (`graftMatches` -> `_querymsf`, lines 387-542)

For each query-family match, constructs a query MSF (Multiple Sequence Format) string that aligns the query into the family's existing MSA using HMM match states.

**`graftMatches` (lines 387-412):**
- Iterates over all matched PANTHER families.
- Skips families without annotations (`$options->{pthrAnn}`).
- Calls `_getAlignLength` to determine MSA width, then `_graftPipeline` per query.

**`_querymsf` (lines 477-542):**
- Initializes a gap-filled array of length `$length_msf`.
- Processes domains sorted by score (highest first).
- For each domain, marks HMM positions as used and fills in the query residues at the corresponding MSA positions.
- Overlapping domains (by HMM position) are skipped with a warning.
- HMM insert states (`.` in the HMM alignment) are skipped -- only match/delete states map to MSA columns.

**`_getAlignLength` (lines 425-455):**
- Reads the first sequence from the `.AN.fasta` file.
- Assumes alignment is wrapped at exactly 80 characters per line.
- Calculates: `(num_lines - 1) * 80 + length(last_line)`.

**Discussion:**
- The alignment length calculation is brittle. It hardcodes 80-character line wrapping and will produce wrong results if the FASTA file uses a different wrap width, or if the last line coincidentally has 80 characters.
- `_querymsf` uses `eq` (string equality) instead of `==` (numeric equality) for length comparison (line 535). This works in Perl but is semantically misleading.
- The overlap detection marks individual HMM positions as used (`@used_positions`), which correctly handles partial overlaps. However, no merging is attempted -- if two high-scoring domains partially overlap, the lower-scoring one is entirely discarded rather than keeping its non-overlapping region.
- The `score` field needed for sorting is only populated by the hmmscan parser (see Stage 2 concern).

### Stage 4: RAxML Grafting (`_generateFasta`, `_runRAxMLAndAnnotate`, lines 548-619)

Writes the query MSF and reference MSA into a combined FASTA file, then runs RAxML to place the query into the reference tree.

**`_generateFasta` (lines 548-570):**
- Writes query sequence (wrapped at 80 chars) with prefix `query_`.
- Appends the entire family `.AN.fasta` alignment file.
- Non-word characters in query ID are replaced with underscores.

**`_runRAxMLAndAnnotate` (lines 574-619):**
- Creates a per-query, per-family temporary directory: `<tmpDir>/<pthr>_<queryid>_raxml<pid>`.
- Runs RAxML in placement mode: `raxmlHPC-SSE3 -f y -p 12345 -m PROTGAMMAWAG -T 4 -G 0.05`.
- Parses the jplace (JSON) output to determine graft placement.
- Falls back to "root" if RAxML fails or produces no placement.

**Discussion:**
- **Fixed random seed** (`-p 12345`) ensures reproducibility but users cannot override it.
- **Hardcoded 4 threads** (`-T 4`) for RAxML regardless of available cores or user's `-cpus` flag (which only applies to HMMER).
- The temp directory cleanup in `_graftPipeline` (line 468) uses `rm -rf $tmpDir/*` via `system()`. This is dangerous: if `$tmpDir` is unset or set to an unexpected value, this could delete unintended files. A Perl-native solution would be safer (the code has a TODO acknowledging this).
- RAxML errors are caught with `try/catch` but only warned about -- the pipeline continues and may produce incorrect "root" annotations for failed placements. There is no way to distinguish between a legitimate root placement and a failure fallback in the output.
- The `queryid` sanitization (`s/[^\w]/\_/g`) is done in both `_generateFasta` and `_runRAxMLAndAnnotate` independently. If either is called without the other, the ID could be inconsistent.

### Stage 5: Annotation Lookup (`mapto`, `commonancestor`, `printResults`, lines 629-806)

Maps RAxML placement to the original PANTHER tree nodes and retrieves annotations.

**`mapto` (lines 670-806):**
- Reads the RAxML `.jplace` JSON output.
- Parses the embedded tree string by extracting `AN####` node labels and their branch indices.
- For single placements: returns the corresponding `AN` node (leaf) or all descendants (internal node).
- For multiple placements: finds the common ancestor of all placement nodes in the RAxML tree, returns its descendants.

**`commonancestor` (lines 631-667):**
- Given a semicolon-separated list of AN nodes and a PANTHER family, finds their lowest common ancestor in the original Newick tree.
- Uses Bio::TreeIO to parse the Newick tree and traverses ancestors.

**Discussion:**
- The `mapto` function is the most complex and fragile part of the codebase. It performs regex surgery on the RAxML tree string, stripping and remapping node labels and branch indices. Any change in RAxML's jplace format would break this.
- The tree string manipulation (`$treestring =~ s/:[0-9\.]+\{([0-9]+)\}/R$1/g`) assumes a specific format from RAxML's jplace output. This format is not formally documented as stable.
- For multiple placements, the common ancestor search uses `@ancestororder` built only from the **last** placement (`$indicator = 1` only when `$i == $nloc-1`). This means the ancestor ordering depends on the last placement's position in the tree, which may not always find the deepest (most specific) common ancestor.
- `mapto` references `Dumper` (line 788) from Data::Dumper but this module is never imported. This will crash on the error path.
- Several `print` and `print STDERR` messages reference "longid" which appears to be a placeholder from an earlier version of the code, not the actual query ID.

## Supporting Scripts

### precompute.pl

Preprocesses raw PANTHER library files (annotation.dat, node.dat, gene_node.dat, annotation_qualifier.dat) into the annotation lookup file (`PAINT_Annotatations_TOTAL.txt`) used by the main pipeline.

**Concerns:**
- References `HAIMING::ToNewick` module (line 8) via hardcoded path (`use lib q(HAIMING)`). This module is not included in the repository.
- Has `#use strict;` commented out (line 10), meaning variable declaration errors will silently pass.
- Uses `@pc` as an array (line 96-99) without `my` declaration -- relies on implicit global, which conflicts with `use warnings`.
- Contains numerous hardcoded paths in comments from previous PANTHER versions.
- Closes `LEAF` and `INTER` filehandles twice (lines 143-144 and 180-181).

### cgi-bin/treeGrafter.cgi

Web service wrapper that accepts a protein sequence via GET/POST, runs the pipeline, and returns XML.

**Concerns:**
- **Command injection vulnerability**: The UUID is used in shell commands (line 127) and while UUIDs themselves are safe, the overall pattern of constructing shell commands from variable interpolation is risky.
- Manually parses URL-encoded form data (lines 82-97) instead of using CGI.pm's built-in parameter handling.
- Hardcoded `hmmscan` algorithm (line 127) -- no option for hmmsearch.
- References a hardcoded `$library` path that is overridden by config, but the initial assignment creates confusion.
- The `usage()` sub describes the program as "myprog.pl" (line 274) -- clearly a template that was never customized.
- Only returns the first result (`last;` on line 148), even if the query matches multiple families.

### sto_to_fasta.py

Converts Stockholm alignment files (`.sto`) to FASTA format (`.AN.fasta`), mapping sequence long IDs to node AN numbers.

**Concerns:**
- Hardcoded dependency on SLURM environment variable `SLURM_PROCID` (line 76) -- cannot run outside SLURM.
- Bare `except:` clause (line 53) catches all exceptions including KeyboardInterrupt.
- Line wrapping is hardcoded to 80 characters (line 107), matching the assumption in `_getAlignLength`.

### gettreesandmsf_hpc.pl

HPC batch script that extracts trees and runs `hmmalign` for PANTHER families.

**Concerns:**
- Uses backtick execution with hardcoded binary path: `` `/home/pmd-02/.../hmmalign -o $sto $hmm $fasta` `` (line 60).
- Depends on `FamLibBuilder` module which is not included in the repository.
- SLURM-dependent (`$ENV{SLURM_PROCID}`).

## Control Flow Summary

```
main
  |-- processOptions()          Parse CLI args, validate paths, load annotations
  |     |-- autodetect()        Choose hmmscan vs hmmsearch (if -auto)
  |
  |-- runhmmer()                Run HMMER (unless precomputed)
  |
  |-- parsehmmer()              Parse HMMER output into $matches hash
  |
  |-- graftMatches()            For each family match:
  |     |-- _getAlignLength()     Read MSA width from .AN.fasta
  |     |-- _graftPipeline()      Per query:
  |           |-- _querymsf()       Build query MSF from HMM alignment
  |           |-- _generateFasta()  Write combined FASTA (query + family MSA)
  |           |-- _runRAxMLAndAnnotate()
  |                 |-- system("raxmlHPC-SSE3 ...")   RAxML placement
  |                 |-- mapto()       Parse jplace, find placement nodes
  |                 |-- commonancestor()  LCA of multiple placement nodes
  |
  |-- printResults()            Write tab-separated output
  |-- rmdir / cleanup
```

## Key Data Structures

### `$options` (hash ref)
Global configuration propagated through all functions. Contains:
- File paths: `fastafile`, `outfile`, `pantherhmm`, `pantherdir`, `tmpDir`, `hmmerout`
- Algorithm: `algo`, `hmmerprecal`, `cpus`, `keep`
- Loaded data: `annotations` (node->annotation map), `pthrAnn` (families with annotations)

### `$matches` (hash ref)
```
$matches->{$pthr_family}->{$query_id} = {
    hmmstart   => [@starts],     # HMM start positions per domain
    hmmend     => [@ends],       # HMM end positions per domain
    score      => [@scores],     # Domain scores (hmmscan only!)
    hmmalign   => [@hmm_seqs],   # HMM consensus alignment strings
    matchalign => [@query_seqs], # Query alignment strings
}
```

### `$allResults` (array ref)
Flat list of tab-separated result strings, one per query.

## Architectural Concerns Summary

1. **Monolithic design**: All pipeline logic in one file with no separation of concerns. Parsing, computation, I/O, and annotation lookup are interleaved. This makes the code difficult to test, extend, or maintain.

2. **No error propagation**: Failures in HMMER, RAxML, or tree parsing are caught and warned about, but the pipeline continues. The output contains no indication of which results are degraded by failures (e.g., "root" fallback).

3. **hmmscan/hmmsearch asymmetry**: The two parsing paths have diverged -- hmmsearch does not extract domain scores, which breaks the multi-domain overlap resolution in `_querymsf`. This is a latent bug.

4. **Shell injection surface**: Multiple `system()` calls construct commands from user-controllable values (file paths, query IDs) without escaping. While file paths are validated for existence, specially crafted filenames could potentially inject shell commands.

5. **Fragile output parsing**: Both HMMER and RAxML outputs are parsed via regex on human-readable output formats rather than machine-readable formats (e.g., HMMER's `--domtblout`, RAxML's structured JSON). This creates tight coupling to specific tool versions.

6. **Missing dependency**: `Data::Dumper` is used in `mapto` error path but never imported. `precompute.pl` requires `HAIMING::ToNewick` which is not in the repo.

7. **Temp file management**: Uses `rm -rf` via `system()` for cleanup, which is dangerous. The tmpDir is shared across pipeline iterations, and `rmdir` at the end (line 68) will fail if any files remain (e.g., from a -k run that was interrupted).

8. **No parallelism within the pipeline**: Each query is processed sequentially through MSA construction and RAxML placement, even though these are independent per query. For large input files, this is a significant bottleneck.

9. **Hardcoded biological model assumptions**: WAG+GAMMA substitution model (`PROTGAMMAWAG`), RAxML random seed, thread count, and the 40x heuristic are all baked in with no user override.

10. **Pre-compiled binary (`treeGrafter101`)**: A macOS universal binary (x86_64 + i386) is checked into the repo. It is 12.9MB, will not run on Linux/ARM, and its provenance is unclear. It predates the current source code state.
