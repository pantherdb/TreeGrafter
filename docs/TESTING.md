# TreeGrafter Testing

## Current Testing Strategy

TreeGrafter has no test framework, no unit tests, and no CI/CD. The only verification method is a manual end-to-end diff:

```bash
perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/test.out -d ./Test/PANTHER_mini -auto
diff ./Test/sample.output ./Test/test.out
```

This compares the pipeline's output against a reference file (`sample.output`) for a fixed input (`sample.fasta`) using a minimal PANTHER dataset (`PANTHER_mini`).

## Test Data Inventory

| File | Size | Purpose |
|---|---|---|
| `Test/sample.fasta` | 36 KB | Input: 125 protein sequences from diverse organisms |
| `Test/sample.output` | 16 KB | Expected output: 125 annotated results |
| `Test/sample.fasta.hmmscan.out` | 1.9 MB | Precomputed HMMER hmmscan results |
| `Test/PANTHER_mini/` | (not in repo) | Minimal PANTHER data directory |

### Key Observation: PANTHER_mini Not in Repository

The `PANTHER_mini` test data directory is **not included in the git repository**. This means:
- The test cannot be run from a fresh clone without additional setup.
- There is no documentation on how to obtain or create `PANTHER_mini`.
- The reference output (`sample.output`) was generated against a specific version of `PANTHER_mini` that may no longer be available.

The precomputed HMMER output file (`sample.fasta.hmmscan.out`) allows partial testing by skipping the HMMER stage:
```bash
perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/test.out -d ./Test/PANTHER_mini \
    -algo hmmscan -hmmer ./Test/sample.fasta.hmmscan.out
```

But this still requires `PANTHER_mini` for the MSA, tree, and annotation files.

## What the Existing Test Covers

The single end-to-end test covers:
- hmmscan output parsing (via precomputed results)
- OR autodetect + HMMER execution (if running with `-auto`)
- MSA construction for all 125 queries
- RAxML placement for queries matching annotated families
- Annotation lookup and output formatting

All 125 queries in the test data match a single family (`PTHR10000`), providing no coverage of:
- Multi-family scenarios
- Queries with no match
- Queries matching families without annotations

## Test Gaps

### No Unit Tests

Every function in `treeGrafter.pl` could benefit from isolated testing:

| Function | Testable Aspects | Risk if Untested |
|---|---|---|
| `autodetect()` | Threshold behavior, edge cases (empty files, single-model HMMs) | Wrong algorithm selected, wasted compute |
| `parsehmmer()` | Both hmmscan and hmmsearch paths, multi-domain handling, edge cases | Silent data loss, wrong family assignment |
| `_querymsf()` | Overlap resolution, insert state handling, length validation | Corrupted MSA, RAxML failure |
| `_getAlignLength()` | Non-80 wrap widths, single-line sequences, empty files | Wrong MSA length, downstream failures |
| `mapto()` | Single/multiple placements, jplace parsing, edge cases | Wrong graft point, wrong annotations |
| `commonancestor()` | Two-node LCA, deep trees, root-level LCA | Overly general annotations |
| `processOptions()` | Missing files, invalid algo, both algo+auto specified | Crash or silent misconfiguration |

### No Coverage of hmmsearch Path

The precomputed test data is hmmscan output. The hmmsearch parsing path (`parsehmmer` else branch) has:
- No test data
- No reference output
- A known bug: does not extract domain scores (see ARCHITECTURE.md)

### No Edge Case Coverage

Untested scenarios include:

1. **Empty input FASTA** -- What happens? (Probably creates empty output, but HMMER may error)
2. **Single sequence input** -- Autodetect will always pick hmmscan; is this correct?
3. **Very long sequences** -- Do HMMER memory limits or RAxML timeout affect results?
4. **Sequences with non-standard amino acids** (B, J, X, Z, *) -- HMMER handles these but does RAxML?
5. **Gene fusion events** -- Query matching two PANTHER families (only first match is used)
6. **Multi-domain queries** -- Overlap resolution logic is untested
7. **Families without annotations** -- Warns to stderr but no output; is this documented behavior?
8. **Missing tree/MSA files** -- Returns 0 from `_getAlignLength`, skips with warning
9. **RAxML failures** -- Falls back to "root" placement silently
10. **Concurrent execution** -- Two pipeline runs on the same input create conflicting temp files
11. **Unicode/special characters in sequence IDs** -- Sanitization to `[^\w]` may behave differently across Perl versions
12. **Extremely large MSAs** -- Memory usage of `_querymsf` (array of individual characters)

### No Regression Testing for Bug Fixes

The commit history shows bug fixes (e.g., "Handle multi-domain overlaps for #7") but there are no test cases that specifically exercise the fixed behavior. Without regression tests, these bugs could easily be reintroduced.

### No Performance Testing

- No benchmarks for scaling behavior (100 vs 10,000 vs 1,000,000 sequences)
- No memory profiling (the entire annotation file is loaded into memory)
- No timing data for the RAxML bottleneck (sequential per-query execution)

## Test Infrastructure Concerns

### Nondeterministic Output

The output order depends on Perl hash iteration order, which is randomized in Perl 5.18+. This means `diff` against `sample.output` may fail even when results are correct, simply because lines are in a different order. A proper comparison would need to sort both files first:

```bash
sort ./Test/sample.output > /tmp/expected_sorted
sort ./Test/test.out > /tmp/actual_sorted
diff /tmp/expected_sorted /tmp/actual_sorted
```

**This is likely already a problem** -- anyone running the test on Perl 5.18+ and getting diff failures may incorrectly conclude the code is broken.

### External Tool Dependencies

Running the full test requires:
- `raxmlHPC-SSE3` (specific binary name, SSE3 instruction set)
- `hmmscan` and/or `hmmsearch` (HMMER 3.1b2)
- Perl with Bio::TreeIO, Try::Tiny, JSON::Parse, IO::String

The Dockerfile pins HMMER 3.1b2 and builds RAxML from source, but local testing has no version pinning. Different HMMER or RAxML versions could produce slightly different placements, causing legitimate diff failures.

### No Test Harness for Precompute Scripts

`precompute.pl`, `sto_to_fasta.py`, `gettreesandmsf_hpc.pl`, and `fix_bifurcate_files.py` have no tests at all. These scripts generate the data that the main pipeline depends on. Bugs in preprocessing would silently corrupt the data directory.

## Recommendations

### Short-Term (Low Effort)

1. **Add sorted diff to test instructions**:
   ```bash
   diff <(sort ./Test/sample.output) <(sort ./Test/test.out)
   ```

2. **Include PANTHER_mini in the repository** (or provide download instructions/script).

3. **Add hmmsearch test data**: Generate `sample.fasta.hmmsearch.out` and a corresponding reference output.

4. **Add a shell-based test runner**:
   ```bash
   #!/bin/bash
   set -e
   perl treeGrafter.pl -f ./Test/sample.fasta -o /tmp/test_hmmscan.out \
       -d ./Test/PANTHER_mini -algo hmmscan -hmmer ./Test/sample.fasta.hmmscan.out
   diff <(sort ./Test/sample.output) <(sort /tmp/test_hmmscan.out) && echo "PASS" || echo "FAIL"
   ```

### Medium-Term (Moderate Effort)

4. **Extract parseable functions for unit testing**: Refactor `parsehmmer`, `_querymsf`, `_getAlignLength`, `mapto`, and `commonancestor` into a Perl module that can be tested with `Test::More`.

5. **Add edge case test data**:
   - Single sequence input
   - Empty input
   - Sequence matching a family without annotations
   - Multi-domain sequence with overlapping domains

6. **Add validation assertions**: At key pipeline stages, validate data structure invariants (e.g., all matched PANTHER families exist in the data directory, MSA lengths match, etc.)

### Long-Term (Significant Effort)

7. **Adopt a test framework**: `Test::More` (Perl) or `pytest` (if migrating to Python) with structured test fixtures.

8. **CI/CD pipeline**: Docker-based CI that runs the full end-to-end test on every commit. The Dockerfile already provides the needed environment.

9. **Fuzzing/property-based testing**: For the parsing functions, use randomized HMMER output to find crash-inducing inputs.

10. **Performance benchmarks**: Track execution time and memory for reference inputs across commits to catch performance regressions.

## Reproducibility Concerns

The test data (`sample.output`) was generated at a specific point in time with specific tool versions. Without recording:
- HMMER version used
- RAxML version used
- Perl version used
- PANTHER data version used
- Operating system

...it is impossible to determine whether a diff failure is a real regression or a tool version artifact. The Dockerfile partially addresses this by pinning HMMER 3.1b2 and building a specific RAxML version, but the PANTHER data version is not pinned.
