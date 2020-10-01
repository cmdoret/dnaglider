package pkg

import (
	"log"
	"math"
	"os"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/unikmer"
)

// KmerProfile stores kmers and their frequencies for a given kmer length
type KmerProfile struct {
	K       int
	Profile map[string]float64
}

// GetSeqKmers compute the k-mer profile of a sequence and increments counts
// in the KmerProfile accordingly
func (p *KmerProfile) GetSeqKmers(seq *seq.Seq) {
	var iter *unikmer.Iterator
	var err error
	var ok bool
	var code uint64

	// Count canonical k-mers for a linear sequence
	iter, err = unikmer.NewKmerIterator(seq, p.K, true, false)
	for {
		code, ok, err = iter.NextKmer()
		if err != nil {
			log.Fatal(err)
			os.Exit(-1)
		}
		if !ok {
			break
		}
		p.Profile[string(unikmer.Decode(code, p.K))]++
	}
}

// CountsToFreqs transforms counts in a KmerProfile into frequencies.
func (p *KmerProfile) CountsToFreqs() {
	var nKmers float64
	// Normalize counts to frequencies
	nKmers = float64(len(p.Profile))
	for kmer, count := range p.Profile {
		p.Profile[kmer] = count / nKmers
	}
}

// KmerDist computes the euclidean distance between a reference k-mer profile
// and another profile. The reference is assumed to include all k-mers
// present in the profile.
func (p *KmerProfile) KmerDist(ref KmerProfile, profile map[string]float64) float64 {
	var dist float64
	for kmer, freq := range ref.Profile {
		dist += math.Pow(freq-profile[kmer], 2.0)
	}
	dist = math.Sqrt(dist)
	return dist
}
