package cli

import (
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"strconv"
	"strings"

	"github.com/cmdoret/dnaglider/dnaglider/pkg"
)

// Command-line flags.
var (
	fasta   = flag.String("fasta", "-", "Input genome. '-' reads from stdin.")
	out     = flag.String("out", "-", "Path to output file. '-' writes to stdout.")
	threads = flag.Int("threads", 1, "Number of CPU threads")
	fields  = flag.String(
		"fields",
		"GC",
		"Statistics to report in output fields. Multiple comma-separated "+
			"values can be provided.\nValid fields are: \n"+
			"\tGC: GC content (0 to 1)\n"+
			"\tGCSKEW: G/C skew (-1 to 1)\n"+
			"\tATSKEW: A/T skew (-1 to 1)\n"+
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
	winSize   = flag.Int("window", 100, "Size of the sliding window.")
	winStride = flag.Int("stride", 100, "Step between windows.")
	version   = flag.Bool("version", false, "Version")
)

// Given a string of comma-separated field names, and a string of
// comma-separated k-mer lengths, build a list of field names. If
// KMER is amon field names, a new field will be generated for each
// k-mer length (e.g. 4MER, 5MER). Otherwise, they will be ignored.
func parseFields(fieldString string, kmerString string) ([]string, []int) {
	metrics := strings.Split(fieldString, ",")
	kmers := strings.Split(kmerString, ",")
	outFields := make([]string, 0)
	outLengths := make([]int, 0)
	for _, metric := range metrics {
		// If k-mers are requested, parse lengths
		if metric == "KMER" {
			for _, k := range kmers {
				kInt, err := strconv.Atoi(k)
				if err != nil {
					log.Fatal(err)
					os.Exit(-1)
				}
				outLengths = append(outLengths, kInt)
				outFields = append(outFields, k+"MER")
			}
		} else {
			outFields = append(outFields, metric)
		}
	}
	return outFields, outLengths
}

// Run reads input genome, chunk it, compute statistics in sliding
// windows and report them. It takes the software version as input.
func Run(semVer string) (err error) {
	var outf io.Writer
	var chunkSize int
	// We'll store the reference profile for each k-mer length
	var kmerLengths []int
	flag.Parse()
	// pprof.StartCPUProfile(os.Stderr)
	// defer pprof.StopCPUProfile()
	// trace.Start(os.Stderr)
	// defer trace.Stop()
	if *version == true {
		fmt.Println(semVer)
		os.Exit(0)
	}
	metrics, kmerLengths := parseFields(*fields, *kmers)
	runtime.GOMAXPROCS(*threads)
	// For each requested kmer length, get the reference profile
	refProfile := make(map[int]pkg.KmerProfile, len(kmerLengths))
	for _, k := range kmerLengths {
		refProfile[k] = pkg.FastaToKmers(*fasta, k)
	}
	chunkSize = 1000
	genome := pkg.StreamGenome(*fasta, 2)
	chunks := pkg.ChunkGenome(genome, *winSize, *winStride, chunkSize)
	results := pkg.ConsumeChunks(chunks, metrics, refProfile)
	// Format each chunk's results into CSV and sends it to an io.writer.
	if *out == "-" {
		outf = os.Stdout
	} else {
		outf, err = os.Create(*out)
		if err != nil {
			log.Fatal("Error opening output file: ", err)
		}
	}
	w := csv.NewWriter(outf)
	w.Comma = '\t'
	res := <-results
	w.Write(res.Header)
	w.WriteAll(res.Data)
	for res := range results {
		w.WriteAll(res.Data)
		if err := w.Error(); err != nil {
			log.Fatal("Error writing csv: ", err)
		}
	}
	return nil
}
