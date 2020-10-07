package pkg

import (
	"io"

	"github.com/shenwei356/bio/seqio/fastx"
)

// StreamGenome reads records one by one from an input fasta file and sends them
// to a channel for downstream processing
func StreamGenome(fasta string, bufSize int) <-chan fastx.Record {
	var record *fastx.Record
	recordChan := make(chan fastx.Record, bufSize)
	reader, err := fastx.NewDefaultReader(fasta)
	// Can't read input path
	checkError(err)
	// Read records asynchronously and send them to a channel
	go func() {
		for {
			record, err = reader.Read()
			if err != nil {
				// Reached last record
				if err == io.EOF {
					close(recordChan)
					return
				}
				checkError(err)
				break
			}
			recordChan <- *record.Clone()
		}
		close(recordChan)
		return
	}()
	return recordChan

}

// FastaToKmers reads all records in a fasta file and computes its k-mer profile
func FastaToKmers(fasta string, k int) KmerProfile {
	var record *fastx.Record
	var profile = KmerProfile{k, make(map[uint64]float64)}
	reader, err := fastx.NewDefaultReader(fasta)
	checkError(err)
	for {
		record, err = reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			checkError(err)
			break
		}
		// Add K-mer counts for each record
		profile.GetSeqKmers(record.Seq)
	}
	// Once all records have been computed, convert to freqs
	profile.CountsToFreqs()
	return profile
}
