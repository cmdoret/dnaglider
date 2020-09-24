package cli

import (
	"flag"

	"github.com/cmdoret/dnaslider/pkg/io"
)

// Command-line flags.
var (
	fasta   = flag.String("fasta", "Path to input genome")
	version = flag.String("version", "0.0.0", "Version")
)

func Run() (err error) {

	flag.Parse()
	genome := io.ReadFasta(*fasta)
}
