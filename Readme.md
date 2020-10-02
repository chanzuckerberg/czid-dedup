# idseq-dedup

[![codecov](https://codecov.io/gh/chanzuckerberg/idseq-dedup/branch/main/graph/badge.svg?token=LMcriTjfuH)](coverage)

idseq-dedup reads FASTA or FASTQ files and outputs versions of those files with duplicate reads removed. A duplicate read in this case is a read that is either identical to another read or shares a prefix of length `-l` with another read. It also outputs a cluster file that lists the representative read ID of each read ID. The representative read ID is the read ID that makes it into the output file, so if a read is found to be a duplicate with a previous read, it will be filtered out of the FASTA/FASTQ output and paired with the read ID of the previous duplicate read in the cluster output file. Representative read IDs are paired with themselves. idseq-dedup can also process paired-end reads. Paired reads are only considered identical if both reads are duplicates to both reads in a previous pair. The order of the input files is preserved. The representative read will always be the first read of its type.

idseq-dedup is still in early development. Binary releases coming soon!
