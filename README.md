# skimmer

Streaming FASTQ quality filter for long reads. Reads four-line FASTQ on stdin (or
from a file), trims a low-quality 3' tail, drops records below a length or mean
Phred threshold, and writes surviving records to stdout.

Built because the filters I had either buffered the whole file, rounded quality
scores in ways that changed downstream basecalling decisions, or pulled in a
dependency tree large enough to make a container rebuild a chore. `skimmer` holds
one record in memory at a time and has no crate dependencies.

## What it does

Per record, in order:

1. Parse and validate: header `@`, separator `+`, sequence length == quality
   length. Malformed records abort with a message naming the line.
2. Trim the 3' tail: from the end, drop every window whose mean Phred score is
   below `--min-mean`. Windows are `--window` bases wide (default 20).
3. Drop the record if the kept sequence is shorter than `--min-len`, or if its
   mean Phred score is below `--min-mean`.
4. Otherwise write the record unchanged apart from the tail trim.

Counters go to stderr; data goes to stdout, so `skimmer ... | gzip > out.fastq.gz`
works without flags.

## Input format

Plain, uncompressed FASTQ only. One record is four lines:

```
@read-name
ACGTACGT
+
IIII#III
```

Quality is interpreted as Phred+33 (Sanger/Illumina 1.8+). CRLF line endings are
accepted. gzip is out of scope; pipe through `gzip -dc` first.

## Install

```
cargo build --release
./target/release/skimmer --help
```

Requires a stable Rust toolchain (edition 2021). No external crates.

## Usage

```
# keep reads of at least 1000 bp whose mean quality clears Q12
gzip -dc longreads.fastq.gz \
  | skimmer --min-len 1000 --min-mean 12 --window 30 \
  | gzip > filtered.fastq.gz

# score only: report how much the file would lose, keep everything but the
# low-quality 3' tails
skimmer --min-mean 10 < in.fastq > trimmed.fastq
```

| flag | default | meaning |
| --- | --- | --- |
| `--min-len` | `0` | drop records shorter than this after trimming |
| `--min-mean` | `0` | drop records below this mean Phred score; `0` disables tail trimming |
| `--window` | `20` | window width used for tail trimming |
| `-h`, `--help` | | print usage |

Exit codes: `0` success, `1` I/O or parse error, `2` bad usage.

## Memory bounds

One record plus two 1 MiB `BufReader`/`BufWriter` buffers. Peak RSS stays flat at
roughly 2 MiB regardless of file size or read length. A single record longer than
the input buffer is read into one `Vec`, so the true bound is
`max(record length, 2 MiB)`.

## Layout

```
src/main.rs      argument parsing, streaming loop, counters
src/fastq.rs     record type, validating reader, writer
src/quality.rs   Phred+33 decoding, mean score, tail trim
```

## Testing

```
cargo test
```

Unit tests in `quality.rs` cover empty input, all-low tails and the no-trim case.
Records are compared byte for byte, so a filter change that alters quality
strings shows up immediately.

## Citation

If this ends up in a methods section, cite the repository:

```
Zhao L. skimmer: streaming FASTQ quality filtering for long reads.
https://github.com/j8myl0dsh0/skimmer-seq (2026).
```

Phred scores follow Ewing & Green, 1998, *Genome Research* 8:186-194; the Q33
encoding is the one described in Cock et al., 2010, *Nucleic Acids Research*
38:e171.

## License

MIT, see `LICENSE`.
