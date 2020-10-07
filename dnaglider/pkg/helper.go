package pkg

import (
	"errors"
	"fmt"
	"os"
	"regexp"
)

var fieldRegex = regexp.MustCompile("^(GC|GCSKEW|ATSKEW|ENTRO|[0-9]+MER)$")

func checkError(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func checkFields(fields []string) error {
	var err error
	for _, field := range fields {
		if !fieldRegex.Match([]byte(field)) {
			err = errors.New("Invalid field name")
			return err
		}
	}
	return err
}
