package pkg

import (
	"fmt"
	"math"
	"os"
	"regexp"

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

// SelectFieldStat will run a function to compute a metric on input sequence
// the metric is selected according to the input field name.
func SelectFieldStat(field string, seq *seq.Seq, ref KmerProfile) float64 {
	kmerRegex := regexp.MustCompile("[0-9]+MER")
	switch field {
	case "GC":
		return SeqGC(seq)
	case "SKEW":
		return SeqGCSkew(seq)
	case "ENTRO":
		return SeqEntropy(seq)
	case kmerRegex.Match(field):
		// TODO: Compute Kmer profile and distance
	default:
		fmt.Printf("Error: Invalid metric: %s.\n", field)
		os.Exit(-1)
	}

}
