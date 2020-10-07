package main

import (
	"fmt"
	"os"

	"github.com/cmdoret/dnaglider/dnaglider/cli"
)

// software version string
var VERSION string

func main() {
	if err := cli.Run(VERSION); err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
}
