package cli

import (
	"flag"
	"fmt"
	"runtime"

	"../pkg"
)

// Command-line flags.
var (
	fasta   = flag.String("fasta", "-", "Path to input genome")
	threads = flag.Int("threads", 1, "Number of CPU threads")
	version = flag.String("version", "0.0.0", "Version")
)

// Run reads input genome, chunk it, compute statistics in sliding
// windows and report them
func Run() (err error) {

	flag.Parse()
	var winSize, chunkSize int
	runtime.GOMAXPROCS(*threads)
	*fasta = "./tests/genome.fa"
	winSize = 10
	chunkSize = 100
	genome := pkg.StreamGenome(*fasta, 3)
	chunks := pkg.ChunkGenome(genome, winSize, chunkSize)
	out := pkg.ConsumeChunks(chunks)
	for chunk := range out {
		fmt.Printf("%s", chunk)
	}
	return nil
}
