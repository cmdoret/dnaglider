package pkg

import (
	"math"

	"github.com/shenwei356/bio/seq"
)

// SeqGC computes the fraction of G or C bases in a sequence (GC content).
func SeqGC(seq *seq.Seq) float64 {
	return seq.GC()
}

// SeqGCSkew computes the GC skew of a DNA sequence.
func SeqGCSkew(seq *seq.Seq) float64 {
	c := seq.BaseContent("c")
	g := seq.BaseContent("g")
	return (g - c) / (g + c)
}

// SeqEntropy computes the Shannon entropy (information entropy)
// of a DNA sequence.
func SeqEntropy(seq *seq.Seq) float64 {
	var prob, entro float64
	seqLen := len(seq.Seq)
	for _, letter := range seq.Alphabet.Letters() {
		prob = seq.BaseContent(string(letter)) / float64(seqLen)
		entro -= prob * math.Log2(prob)
	}
	return entro
}

// TODO: Add a function to compute vectors of k-mer frequencies.
// These k-mer frequencies need to be compared against a reference
// profile and some distance metric can then be computed between each
// window's vector and the reference.
