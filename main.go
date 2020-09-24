package main

import (
	"flag"
	"fmt"

	"github.com/cmdoret/dnaslider/cli"
)

// Command-line flags.
var (
	fasta   = flag.String("fasta", "Path to input genome")
	version = flag.String("version", "0.0.0", "Version")
)

func main() {
	if err := cli.Run(); err != nil {
		fmt.Println(err)
	}
}
