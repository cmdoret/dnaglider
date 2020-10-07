package pkg

import (
	"math"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/unikmer"
)

// KmerProfile stores kmers and their frequencies for a given kmer length
type KmerProfile struct {
	K       int
	Profile map[string]float64
}

// NewKmerProfile is a helper function to generate an empty Kmer profile
func NewKmerProfile(k int) KmerProfile {
	return KmerProfile{k, make(map[string]float64)}
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
	checkError(err)
	for {
		code, ok, err = iter.NextKmer()
		checkError(err)
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
func (p *KmerProfile) KmerDist(ref KmerProfile) float64 {
	var dist float64
	for kmer, freq := range ref.Profile {
		dist += math.Pow(freq-p.Profile[kmer], 2.0)
	}
	dist = math.Sqrt(dist)
	return dist
}
