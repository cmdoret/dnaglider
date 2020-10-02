package pkg

import (
	"fmt"
	"math"
	"os"
	"regexp"
	"strconv"

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
	for _, letter := range []string{"A", "C", "T", "G"} {
		prob = seq.BaseContent(letter)
		if prob > 0 {
			entro += (prob * math.Log2(prob))
		}
	}
	return -entro
}

// SelectFieldStat will run a function to compute a metric on input sequence
// the metric is selected according to the input field name.
func SelectFieldStat(field string, seq *seq.Seq, ref map[int]KmerProfile) (float64, error) {
	var err error
	kmerRegex := regexp.MustCompile("([0-9]+)MER")
	switch field {
	case "GC":
		return SeqGC(seq), err
	case "SKEW":
		return SeqGCSkew(seq), err
	case "ENTRO":
		return SeqEntropy(seq), err
	default:
		// Get the k-mer length from field name, compute the
		// sequence k-mer profile and its distance to the ref profile
		if kmerRegex.Match([]byte(field)) {
			k, err := strconv.Atoi(kmerRegex.FindStringSubmatch(field)[1])
			prof := KmerProfile{k, make(map[string]float64)}
			return prof.KmerDist(ref[k]), err
		}
		fmt.Printf("Error: Invalid metric: %s.\n", field)
		os.Exit(-1)
	}
	return 0.0, err
}
