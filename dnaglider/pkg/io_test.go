package pkg

import (
	"testing"
)

var TESTFILE string = "../../tests/genome.fa"

func TestStreamGenome(t *testing.T) {
	var nRecords int
	records := StreamGenome(TESTFILE, 3)
	for rec := range records {
		nRecords++
		if len(rec.Seq.Seq) != 1200 {
			t.Errorf("StreamGenome found incorrect sequence length.")
		}

	}
	if nRecords != 16 {
		t.Errorf("StreamGenome read %d records instead of %d", nRecords, 16)

	}
}
