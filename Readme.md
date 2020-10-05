# idseq-dedup

![Rust](https://github.com/chanzuckerberg/idseq-dedup/workflows/CICD/badge.svg) [![codecov](https://codecov.io/gh/chanzuckerberg/idseq-dedup/branch/main/graph/badge.svg?token=LMcriTjfuH)](coverage) [![GitHub license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://github.com/chanzuckerberg/idseq-web/blob/master/LICENSE) ![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)

idseq-dedup reads FASTA or FASTQ files and outputs versions of those files with duplicate reads removed. A duplicate read in this case is a read that is either identical to another read or shares a prefix of length `-l` with another read. It also outputs a cluster file that lists the representative read ID of each read ID. The representative read ID is the read ID that makes it into the output file, so if a read is found to be a duplicate with a previous read, it will be filtered out of the FASTA/FASTQ output and paired with the read ID of the previous duplicate read in the cluster output file. Representative read IDs are paired with themselves. idseq-dedup can also process paired-end reads. Paired reads are only considered identical if both reads are duplicates to both reads in a previous pair. The order of the input files is preserved. The representative read will always be the first read of its type.

FASTA/FASTQ parsing provided by [rust-bio](https://github.com/rust-bio/rust-bio).

## Installation

### Binary

We release binaries for Linux, MacOS, and Windows. To install one, download the appropriate binary for your operating system from one of our [releases](https://github.com/chanzuckerberg/idseq-dedup/releases/).

### From Source

1. [Install rust/cargo](https://www.rust-lang.org/tools/install) if you haven't already
1. `git clone https://github.com/chanzuckerberg/idseq-dedup.git`
1. `cd idseq-dedup`
1. `cargo build --release`
1. Your executable will be at `idseq-dedup/target/release/idseq-dedup` (with `.exe` if you're on windows)

## Usage

Run:

```bash
idseq-dedup --help
```

for usage information:

```
USAGE:
    idseq-dedup [OPTIONS] --deduped-outputs <deduped-outputs>... --inputs <inputs>...

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --cluster-output <cluster-output>         Output cluster file [default: clusters.csv]
    -o, --deduped-outputs <deduped-outputs>...    Output deduped FASTQ
    -i, --inputs <inputs>...                      Input FASTQ
    -l, --prefix-length <prefix-length>           Length of the prefix to consider
```

### Example Usage

Deduplicate a single end FASTA:

```bash
idseq-dedup -i my-fasta.fasta -o my-deduped-fasta.fasta
```

Deduplicate a single end FASTQ (same as FASTA)

```bash
idseq-dedup -i my-fastq.fastq -o my-deduped-fastq.fastq
```

Deduplicate paired reads (note, inputs are paired to outputs by order not name):

```bash
idseq-dedup \
	-i my-fasta-r1.fasta \
	-i my-fasta-r2.fasta \
	-o my-deduped-fasta-r1.fasta \
	-o my-deduped-fasta-r2.fasta
```

Deduplicate only considering a prefix of length `70`:

```bash
idseq-dedup -l 70 -i my-fasta.fasta -o my-deduped-fasta.fasta
```

Custom cluster file name:

```bash
idseq-dedup -i my-fasta.fasta -o my-deduped-fasta.fasta -c custom-cluster.csv
```

