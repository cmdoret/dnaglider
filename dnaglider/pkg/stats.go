package pkg

import (
	"math"

	"github.com/shenwei356/bio/seq"
)

var fieldDispatcher = func() map[string]func(seq *seq.Seq) float64 {
	fd := make(map[string]func(seq *seq.Seq) float64)
	fd["GC"] = SeqGC
	fd["GCSKEW"] = SeqGCSkew
	fd["ATSKEW"] = SeqATSkew
	fd["ENTRO"] = SeqEntropy
	return fd
}()

// SeqGC computes the fraction of G or C bases in a sequence (GC content).
func SeqGC(seq *seq.Seq) float64 {
	return seq.GC()
}

// SeqGCSkew computes the GC skew of a DNA sequence.
func SeqGCSkew(seq *seq.Seq) float64 {
	c := seq.BaseContent("c")
	g := seq.BaseContent("g")
	if g+c == 0 {
		return math.NaN()
	}
	return (g - c) / (g + c)
}

// SeqATSkew computes the AT skew of a DNA sequence.
func SeqATSkew(seq *seq.Seq) float64 {
	a := seq.BaseContent("a")
	t := seq.BaseContent("t")
	if a+t == 0 {
		return math.NaN()
	}
	return (a - t) / (a + t)
}

// SeqEntropy computes the Shannon entropy (information entropy)
// of a DNA sequence.
func SeqEntropy(seq *seq.Seq) float64 {
	var prob, entro float64
	for _, letter := range []string{"A", "C", "T", "G"} {
		prob = seq.BaseContent(letter)
		if prob > 0 {
			entro += (prob * math.Log2(prob))
		}
	}
	return -entro
}

// SeqKmerDiv will compute the Kmer profile of the input profile
// and compute its distance to a reference k-mer profile.
func SeqKmerDiv(seq *seq.Seq, ref KmerProfile) float64 {
	// Get the k-mer length from field name, compute the
	// sequence k-mer profile and its distance to the ref profile
	prof := NewKmerProfile(ref.K)
	prof.GetSeqKmers(seq)
	prof.CountsToFreqs()
	return prof.KmerCosDist(ref)
}
