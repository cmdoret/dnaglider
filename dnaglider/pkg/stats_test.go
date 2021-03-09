package pkg

import (
	"math"
	"testing"

	"github.com/shenwei356/bio/seq"
)

// TODO: Add test cases with ambiguous or invalid bases

func TestSeqGC(t *testing.T) {
	gcExamples := map[string]float64{"AAAAAA": 0.0, "AAAC": 0.25, "CAAG": 0.5, "GCGC": 1.0}
	for exSeq, expGC := range gcExamples {
		tmpSeq, _ := seq.NewSeq(seq.DNA, []byte(exSeq))
		obsGC := SeqGC(tmpSeq)
		if obsGC != expGC {
			t.Errorf("Wrong GC content for %s: got %f instead of %f", exSeq, obsGC, expGC)
		}
	}
}

func TestGCSkew(t *testing.T) {
	skewExamples := map[string]float64{"AAAC": -1, "CAAG": 0, "GCGG": 0.5}
	// Test valid sequences
	for exSeq, expSkew := range skewExamples {
		tmpSeq, _ := seq.NewSeq(seq.DNA, []byte(exSeq))
		obsSkew := SeqGCSkew(tmpSeq)
		if obsSkew != expSkew {
			t.Errorf("Wrong GC skew value for %s: got %f instead of %f", exSeq, obsSkew, expSkew)
		}
	}
	// Test invalid case
	invalidSeq, _ := seq.NewSeq(seq.DNA, []byte("AAAAAA"))
	if !math.IsNaN(SeqGCSkew(invalidSeq)) {
		t.Errorf("GC skew should have been invalid for: %s", invalidSeq)
	}

}

func TestATSkew(t *testing.T) {
	skewExamples := map[string]float64{"CCCT": -1, "ACCT": 0, "ATAA": 0.5}
	// Test valid sequences
	for exSeq, expSkew := range skewExamples {
		tmpSeq, _ := seq.NewSeq(seq.DNA, []byte(exSeq))
		obsSkew := SeqATSkew(tmpSeq)
		if obsSkew != expSkew {
			t.Errorf("Wrong AT skew value for %s: got %f instead of %f", exSeq, obsSkew, expSkew)
		}
	}
	// Test invalid case
	invalidSeq, _ := seq.NewSeq(seq.DNA, []byte("CCCCCC"))
	if !math.IsNaN(SeqATSkew(invalidSeq)) {
		t.Errorf("AT skew should have been invalid for: %s", invalidSeq)
	}
}
func TestSeqEntropy(t *testing.T) {
	entroExamples := map[string]float64{"CCCC": 0, "ACGT": 2, "ATTA": 1}
	// Test valid sequences
	for exSeq, expSkew := range entroExamples {
		tmpSeq, _ := seq.NewSeq(seq.DNA, []byte(exSeq))
		obsSkew := SeqEntropy(tmpSeq)
		if obsSkew != expSkew {
			t.Errorf("Wrong entropy value for %s: got %f instead of %f", exSeq, obsSkew, expSkew)
		}
	}

}

func TestSeqKmerDiv(t *testing.T) {
	refSeq, _ := seq.NewSeq(seq.DNA, []byte("GCGCG"))
	querySeq, _ := seq.NewSeq(seq.DNA, []byte("AATAA"))

	refProf := NewKmerProfile(3)
	refProf.GetSeqKmers(refSeq)
	refProf.CountsToFreqs()
	expDist := 1.0
	obsDist := SeqKmerDiv(querySeq, refProf)
	if expDist != obsDist {
		t.Errorf(
			"Wrong kmer distance with reference: got %f instead of %f",
			obsDist,
			expDist,
		)
	}
}
