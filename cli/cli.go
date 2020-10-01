package cli

import (
	"flag"
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
	winSize = 1000
	chunkSize = 1000
	genome := pkg.StreamGenome(*fasta, 3)
	chunks := pkg.ChunkGenome(genome, winSize, chunkSize)
	pkg.ConsumeChunks(chunks)
	return nil
}
