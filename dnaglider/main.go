package main

import (
	"fmt"

	"github.com/cmdoret/dnaglider/dnaglider/cli"
)

// software version string
var VERSION string

func main() {
	if err := cli.Run(VERSION); err != nil {
		fmt.Println(err)
	}
}
