package main

import (
	"fmt"

	"github.com/cmdoret/dnaglider/dnaglider/cli"
)

func main() {
	if err := cli.Run(); err != nil {
		fmt.Println(err)
	}
}
