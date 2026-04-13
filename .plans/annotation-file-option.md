# Annotation File Option (`-a`) Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Allow users to specify the PAINT annotation file path independently from the `-d` PANTHER data directory.

**Architecture:** Add a new `-a` CLI option that accepts an annotation file path. When provided, it overrides the default `$directory/PAINT_Annotations/PAINT_Annotatations_TOTAL.txt`. When omitted, behavior is unchanged (backwards-compatible).

**Tech Stack:** Perl, Getopt::Long

---

## File Structure

- **Modify:** `treeGrafter.pl` — Add `-a` option to GetOptions, thread through processOptions, use as override for annotation file path.

---

### Task 1: Add the `-a` CLI option and wire it through

**Files:**
- Modify: `treeGrafter.pl:15-26` (variable declarations)
- Modify: `treeGrafter.pl:28-42` (GetOptions block)
- Modify: `treeGrafter.pl:46-48` (processOptions call)
- Modify: `treeGrafter.pl:74-76` (processOptions signature)
- Modify: `treeGrafter.pl:127-128` (annotation file path assignment)
- Modify: `treeGrafter.pl:845-867` (usage text)

- [ ] **Step 1: Add `$annotationfile` variable declaration**

At `treeGrafter.pl:15`, add `$annotationfile` to the `my (...)` declaration list:

```perl
my ($fastafile,
    $outfile,
    $raxmlloc,
    $hmmerloc,
    $directory,
    $tmpDir,
    $keep,
    $algo,
    $cpus,
    $auto,
    $hmmer,
    $annotationfile,
    $help);
```

- [ ] **Step 2: Add `-a` to GetOptions**

At `treeGrafter.pl:28-42`, add a new line to the GetOptions block:

```perl
"a=s" => \$annotationfile, # -a for the PAINT annotation file path
```

- [ ] **Step 3: Pass `$annotationfile` to processOptions**

At `treeGrafter.pl:46-48`, update the processOptions call:

```perl
processOptions( $options, $help, $outfile, $directory, $tmpDir, $fastafile,
               $raxmlloc, $hmmerloc, $algo, $auto, $cpus, $keep, $annotationfile);
```

- [ ] **Step 4: Accept `$annotationfile` in processOptions signature**

At `treeGrafter.pl:74-76`, update the sub signature:

```perl
sub processOptions {
  my ( $options, $help, $outfile, $directory, $tmpDir, $fastafile,
        $raxmlloc, $hmmerloc, $algo, $auto, $cpus, $keep, $annotationfile) = @_;
```

- [ ] **Step 5: Use `$annotationfile` with fallback to default**

At `treeGrafter.pl:127-128`, replace the hard-coded assignment with a conditional:

```perl
  # my $annotationfile = "$directory/PANTHER12_PAINT_Annotations/PANTHER12_PAINT_Annotatations_TOTAL.txt";
  if (!$annotationfile) {
    $annotationfile = "$directory/PAINT_Annotations/PAINT_Annotatations_TOTAL.txt";
  }
```

Note: `$annotationfile` is already declared in the outer scope and passed in as a parameter — this shadows appropriately within processOptions. The variable passed in from GetOptions will be `undef` if `-a` was not provided, so the fallback triggers.

- [ ] **Step 6: Update usage text**

At `treeGrafter.pl:852-867`, add the `-a` option to the usage output:

```perl
    -a for the PAINT annotation file path (default: <-d dir>/PAINT_Annotations/PAINT_Annotatations_TOTAL.txt)
```

Insert after the `-d` line.

- [ ] **Step 7: Test with default behavior (no `-a` flag)**

Run:
```bash
perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/test_default.out -d ./Test/PANTHER_mini -auto
diff ./Test/sample.output ./Test/test_default.out
```
Expected: No diff — output matches reference. Confirms backwards compatibility.

- [ ] **Step 8: Test with explicit `-a` flag pointing to the same file**

Run:
```bash
perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/test_explicit.out -d ./Test/PANTHER_mini -auto -a ./Test/PANTHER_mini/PAINT_Annotations/PAINT_Annotatations_TOTAL.txt
diff ./Test/sample.output ./Test/test_explicit.out
```
Expected: No diff — output matches reference. Confirms `-a` flag works.

- [ ] **Step 9: Test with nonexistent `-a` path**

Run:
```bash
perl treeGrafter.pl -f ./Test/sample.fasta -o ./Test/test_bad.out -d ./Test/PANTHER_mini -auto -a /nonexistent/file.txt 2>&1
```
Expected: Dies with "The PANTHER annotation file, /nonexistent/file.txt, does not exist."

- [ ] **Step 10: Commit**

```bash
git add treeGrafter.pl
git commit -m "feat: add -a option for specifying annotation file path separately from -d"
```

---

## Update Documentation

### Task 2: Update CLAUDE.md and docs

**Files:**
- Modify: `CLAUDE.md`
- Modify: `docs/ARCHITECTURE.md` (if it documents CLI flags)

- [ ] **Step 1: Update CLAUDE.md**

In the "Running" section, add `-a` to the key flags list:

```
Key flags: `-f` input FASTA, `-o` output, `-d` PANTHER data dir, `-a` annotation file (overrides default in `-d`), `-algo hmmscan|hmmsearch`, `-auto` auto-select algorithm, `-hmmer` precomputed results, `-k` keep temp files, `-t` tmpdir, `-cpus` HMMER threads.
```

- [ ] **Step 2: Commit docs**

```bash
git add CLAUDE.md
git commit -m "docs: document -a annotation file option"
```
