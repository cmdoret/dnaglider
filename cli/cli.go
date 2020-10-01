package cli

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"strconv"
	"strings"

	"../pkg"
)

// Command-line flags.
var (
	fasta   = flag.String("fasta", "-", "Input genome. '-' reads from stdin.")
	threads = flag.Int("threads", 1, "Number of CPU threads")
	fields  = flag.String(
		"fields",
		"GC",
		"Statistics to report in output fields. Multiple comma-separated "+
			"values can be provided.\nValid fields are: \n"+
			"\tGC: GC content (0 to 1)\n"+
			"\tSKEW: G/C skew (-1 to 1)\n"+
			"\tENTRO: Information entropy of the sequence (0 to 1)\n"+
			"\tKMER: K-mer divergence from the reference (euclidean distance)\n",
	)
	kmers = flag.String(
		"kmers",
		"4",
		"Report k-mer divergence from the genome for the following k-mer "+
			"lengths. Multiple comma separated values can be provided. This only "+
			"has an effect if KMER is specified in -fields.",
	)
	winSize = flag.Int("window", 100, "Size of the sliding window.")
	version = flag.String("version", "0.0.0", "Version")
)

// Run reads input genome, chunk it, compute statistics in sliding
// windows and report them
func Run() (err error) {
	var winSize, chunkSize int
	// We'll store the reference profile for each k-mer length
	var refProfile map[int]pkg.KmerProfile
	var kmerLengths []int
	flag.Parse()
	runtime.GOMAXPROCS(*threads)
	// Parse fields to compute
	metrics := strings.Split(*fields, ",")
	for _, metric := range metrics {
		// If k-mers are requested, parse lengths
		if metric == "KMER" {
			for _, k := range strings.Split(*kmers, ",") {
				kInt, err := strconv.Atoi(k)
				if err != nil {
					log.Fatal(err)
					os.Exit(-1)
				}
				kmerLengths = append(kmerLengths, kInt)
			}
		}
	}
	// For each requested length, get the reference k-mer profile
	for _, k := range kmerLengths {
		refProfile[k] = pkg.FastaToKmers(*fasta, k)
	}
	*fasta = "./tests/genome.fa"
	chunkSize = 1000
	genome := pkg.StreamGenome(*fasta, 3)
	chunks := pkg.ChunkGenome(genome, winSize, chunkSize)
	out := pkg.ConsumeChunks(chunks, metrics, refProfile)
	for chunk := range out {
		fmt.Printf("%s", chunk)
	}
	return nil
}
