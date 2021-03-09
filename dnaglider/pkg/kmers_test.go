package pkg

import (
	"fmt"
	"math"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/unikmer"
)

func TestNewKmerProfile(t *testing.T) {
	tProf := NewKmerProfile(3)
	if tProf.K != 3 {
		t.Errorf("KmerProfile generated with wrong K value.")
	}
}
func TestGetSeqKmers(t *testing.T) {
	testSeq, _ := seq.NewSeq(seq.DNA, []byte("CCTAAA"))
	obsProf := NewKmerProfile(3)
	obsProf.GetSeqKmers(testSeq)
	for code, freq := range obsProf.Profile {
		fmt.Println(string(unikmer.Decode(code, 3)), freq)
	}
	k1, _ := unikmer.Encode([]byte("AGG"))
	k2, _ := unikmer.Encode([]byte("CTA"))
	expProf := map[uint64]float64{k1: 1.0, k2: 1.0}
	if len(obsProf.Profile) != len(expProf) {
		t.Errorf(
			"Number of kmers in profile incorrect: got %d instead of %d",
			len(obsProf.Profile),
			len(expProf),
		)
	}
	for kmer, count := range expProf {
		if obsProf.Profile[kmer] != count {
			t.Errorf(
				"Wrong count for kmer %s: got %f instead of %f",
				unikmer.Decode(kmer, obsProf.K),
				obsProf.Profile[kmer],
				count,
			)
		}
	}
}
func TestCountsToFreqs(t *testing.T) {
	testSeq, _ := seq.NewSeq(seq.DNA, []byte("CCTAAA"))
	obsProf := NewKmerProfile(3)
	obsProf.GetSeqKmers(testSeq)
	obsProf.CountsToFreqs()
	k1, _ := unikmer.Encode([]byte("AGG"))
	k2, _ := unikmer.Encode([]byte("CTA"))
	expProf := map[uint64]float64{k1: 0.5, k2: 0.5}
	for kmer, count := range expProf {
		if obsProf.Profile[kmer] != count {
			t.Errorf(
				"Wrong frequency for kmer %s: got %f instead of %f",
				unikmer.Decode(kmer, obsProf.K),
				obsProf.Profile[kmer],
				count,
			)
		}
	}

}
func TestKmerEuclDist(t *testing.T) {
	testSeq, _ := seq.NewSeq(seq.DNA, []byte("GCGCGC"))
	testRef, _ := seq.NewSeq(seq.DNA, []byte("CCTAAA"))
	profSeq := NewKmerProfile(3)
	profRef := NewKmerProfile(3)
	profSeq.GetSeqKmers(testSeq)
	profRef.GetSeqKmers(testRef)
	profSeq.CountsToFreqs()
	profRef.CountsToFreqs()
	expDist := math.Sqrt(0.5)
	obsDist := profSeq.KmerEuclDist(profRef)
	if expDist != obsDist {
		t.Errorf(
			"Incorrect Euclidean distance between kmer profiles:  got %f instead of %f",
			obsDist,
			expDist,
		)
	}

}

func TestKmerCosDist(t *testing.T) {
	testSeq, _ := seq.NewSeq(seq.DNA, []byte("CCTAAA"))
	testRef, _ := seq.NewSeq(seq.DNA, []byte("CCTAAA"))
	profSeq := NewKmerProfile(3)
	profRef := NewKmerProfile(3)
	profSeq.GetSeqKmers(testSeq)
	profRef.GetSeqKmers(testRef)
	profSeq.CountsToFreqs()
	profRef.CountsToFreqs()
	expDist := 0.0
	obsDist := profSeq.KmerCosDist(profRef)
	if math.Abs(expDist-obsDist) > 0.0000001 {
		fmt.Println()
		t.Errorf(
			"Incorrect cosine distance between kmer profiles:  got %f instead of %f",
			obsDist,
			expDist,
		)
	}

}
