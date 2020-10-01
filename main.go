package main

import (
	"fmt"

	"./cli"
)

func main() {
	if err := cli.Run(); err != nil {
		fmt.Println(err)
	}
}
