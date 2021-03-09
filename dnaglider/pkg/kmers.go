package pkg

import (
	"math"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/unikmer"
)

// KmerProfile stores kmers and their frequencies for a given kmer length
type KmerProfile struct {
	K       int
	Profile map[uint64]float64
}

// NewKmerProfile is a helper function to generate an empty Kmer profile
func NewKmerProfile(k int) KmerProfile {
	return KmerProfile{k, make(map[uint64]float64, int(math.Pow(4, float64(k))))}
}

// GetSeqKmers compute the k-mer profile of a sequence and increments counts
// in the KmerProfile accordingly
func (p *KmerProfile) GetSeqKmers(sq *seq.Seq) {
	var iter *unikmer.Iterator
	var ok bool
	var code uint64
	// Count canonical k-mers for a linear sequence
	iter, _ = unikmer.NewKmerIterator(sq, p.K, true, false)
	for {
		code, ok, _ = iter.NextKmer()
		if !ok {
			break
		}
		p.Profile[code]++
	}
}

// CountsToFreqs transforms counts in a KmerProfile into frequencies.
func (p *KmerProfile) CountsToFreqs() {
	var nKmers float64
	// Compute total number of k-mers
	for _, count := range p.Profile {
		nKmers += count
	}
	// Normalize counts to frequencies
	for kmer, count := range p.Profile {
		p.Profile[kmer] = count / nKmers
	}
}

// KmerEuclDist computes the euclidean distance between a reference k-mer profile
// and another profile. The reference is assumed to include all k-mers
// present in the profile.
func (p *KmerProfile) KmerEuclDist(ref KmerProfile) float64 {
	var dist float64
	for kmer, freq := range ref.Profile {
		// This works because fetching missing k-mer
		// returns zero value for float

		dist += math.Pow(freq-p.Profile[kmer], 2.0)
	}
	dist = math.Sqrt(dist)
	return dist
}

// KmerCosDist computes the cosine distance between a reference k-mer profile
// and another profile. The reference is assumed to include all k-mers
// present in the profile.
func (p *KmerProfile) KmerCosDist(ref KmerProfile) float64 {
	var dist, numerator, denominator, normL, normR float64
	for kmer, freq := range ref.Profile {
		// This works because fetching missing k-mer
		// returns zero value for float
		numerator += freq * p.Profile[kmer]
		normL += math.Pow(freq, 2.0)
		normR += math.Pow(p.Profile[kmer], 2.0)
	}
	denominator = math.Sqrt(normL) * math.Sqrt(normR)
	// Prevent nans
	if denominator == 0.0 {
		denominator = 1.0
	}
	dist = 1.0 - numerator/denominator
	return dist
}
