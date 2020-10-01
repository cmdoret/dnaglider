package pkg

import (
	"testing"
)

func TestStreamGenome(t *testing.T) {
	var nRecords int
	testFile := "../tests/genome.fa"
	records := StreamGenome(testFile, 3)
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
