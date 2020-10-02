## dnaglider

Command line utility to compute sliding window genome statistics from a fasta file.

### Installation:

If [go is installed](https://golang.org/doc/install) on the machine, the program can be built from source using: 

```bash
go install github.com/cmdoret/dnaglider/dnaglider
```

Otherwise, binaries can be downloaded from the github repository.

### Usage:

dnaglider only requires a genome. You can also select a window size and what metrics to compute. For example to compute GC content and GC skew on 8 threads:

```bash
dnaglider -window 1000 -threads 8 -fields "GC,GCSKEW" -fasta ./mygenome.fasta -out gc_stats.tsv
```

Instead of working with input / output files, the program reads from stdin and write to stdout by deault:

```bash
some command genome.fa | dnaglider -fields "GC,GCSKEW" | grep "chr10" > gc_stats_chr10.tsv
```

```
Usage of dnaglider:
  -fasta string
        Input genome. '-' reads from stdin. (default "-")
  -fields string
        Statistics to report in output fields. Multiple comma-separated values can be provided.
        Valid fields are: 
                GC: GC content (0 to 1)
                GCSKEW: G/C skew (-1 to 1)
                ATSKEW: A/T skew (-1 to 1)
                ENTRO: Information entropy of the sequence (0 to 1)
                KMER: K-mer divergence from the reference (euclidean distance)
         (default "GC")
  -kmers string
        Report k-mer divergence from the genome for the following k-mer lengths. Multiple comma separated values can be provided. This only has an effect if KMER is specified in -fields. (default "4")
  -out string
        Path to output file. '-' writes to stdout. (default "-")
  -threads int
        Number of CPU threads (default 1)
  -version string
        Version (default "0.0.0")
  -window int
        Size of the sliding window. (default 100)
```

### Output:

The output files are tab-separated text files with one row per window. The first 3 columns indicate 1-based genomic coordinates and the following column contain statistics computed on the genome.
