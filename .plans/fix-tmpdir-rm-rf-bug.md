# Fix: tmpDir rm -rf catastrophic deletion bug

## Problem

In the `fullgo-paint-update-77-annot-file-option` branch, the `-t` (tmpDir) command-line option was added but has no default value or required-argument check. When `-t` is not passed:

1. `$tmpDir` is `undef`
2. Perl stringifies `undef` to `""` in the concatenation at `_graftPipeline`:
   ```perl
   my $command = "rm -rf ".$options->{tmpDir}."/*";
   ```
3. This produces `rm -rf /*`
4. Since the container runs as root with bind-mounted host directories, this deletes all files in every mounted host directory

## Affected code

- **treeGrafter.pl line ~464** (`_graftPipeline`): `rm -rf` cleanup of temp files
- **treeGrafter.pl line ~66** (main body): `rmdir($options->{tmpDir})` cleanup
- **treeGrafter.pl line ~52-55** (`processOptions`): `tmpDir` assignment and mkdir

## Fix steps

### Step 1: Add default value for tmpDir in processOptions

In `processOptions`, after the existing directory validation, add a default for `tmpDir` when not specified. Fall back to `$directory/tmp` (the original behavior before the branch diverged):

```perl
if (!defined($tmpDir) || $tmpDir eq '') {
  $tmpDir = "$directory/tmp";
}
$options->{tmpDir} = $tmpDir;
if (!-d "$tmpDir") {
  mkdir("$tmpDir");
}
```

### Step 2: Add guard in _graftPipeline before rm -rf

Add a safety check before the `rm -rf` to ensure `tmpDir` is a non-empty, valid path:

```perl
unless($options->{keep}){
  my $tmpDir = $options->{tmpDir};
  if (defined($tmpDir) && $tmpDir ne '' && $tmpDir ne '/' && -d $tmpDir) {
    my $command = "rm -rf " . $tmpDir . "/*";
    system($command);
  }
}
```

### Step 3: Add same guard around rmdir in main body

```perl
if (!$options->{keep}) {
  my $tmpDir = $options->{tmpDir};
  if (defined($tmpDir) && $tmpDir ne '' && $tmpDir ne '/') {
    rmdir($tmpDir);
  }
}
```

### Step 4: Update Dockerfile ENTRYPOINT or add -t to example CMD

Update the example CMD comment in the Dockerfile to include `-t /tmp`:

```dockerfile
# Example CMD
# docker run --rm --name treegrafter -v /path/to/data:/data treegrafter -f ./Test/sample.fasta -o /tmp/sample.1.out -d /data/PANTHER -t /tmp -auto
```

## Files to modify

1. `treeGrafter.pl` — steps 1-3 (on branch `fullgo-paint-update-77-annot-file-option`)
2. `Dockerfile` — step 4

## Testing

- Run with `-t /tmp` explicitly: verify temp files go to `/tmp` and cleanup works
- Run WITHOUT `-t`: verify it defaults to `$directory/tmp` and does NOT produce `rm -rf /*`
- Run with `-k` flag: verify no cleanup occurs
