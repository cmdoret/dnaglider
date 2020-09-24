package dnaslider

import (
	"io"
	"log"
	"os"

	"github.com/shenwei356/bio/seqio/fastx"
)

// StreamGenome reads records one by one from an input fasta file and sends them
// to a channel for downstream processing
func StreamGenome(Fasta string) <-chan fastx.Record {
	var record *fastx.Record
	recordChan := make(chan fastx.Record)
	reader, err := fastx.NewDefaultReader(Fasta)
	// Can't read input path
	if err != nil {
		log.Fatal(err)
		os.Exit(-1)
	}
	// Read records asynchronously and send them to a channel
	go func() {
		for {
			record, err = reader.Read()
			if err != nil {
				// Reached last record
				if err == io.EOF {
					break
				}
				break
			}
			recordChan <- *record
		}
		close(recordChan)
		return
	}()
	return recordChan

}
