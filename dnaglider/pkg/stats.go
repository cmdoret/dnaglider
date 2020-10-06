package pkg

import (
	"fmt"
	"math"
	"regexp"
	"strconv"

	"github.com/shenwei356/bio/seq"
)

var kmerRegex = regexp.MustCompile("([0-9]+)MER")

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

// SelectFieldStat will run a function to compute a metric on input sequence
// the metric is selected according to the input field name.
func SelectFieldStat(field string, seq *seq.Seq, ref map[int]KmerProfile) (float64, error) {
	var err error
	switch field {
	case "GC":
		return SeqGC(seq), err
	case "GCSKEW":
		return SeqGCSkew(seq), err
	case "ATSKEW":
		return SeqATSkew(seq), err
	case "ENTRO":
		return SeqEntropy(seq), err
	default:
		// Get the k-mer length from field name, compute the
		// sequence k-mer profile and its distance to the ref profile
		if kmerRegex.Match([]byte(field)) {
			k, err := strconv.Atoi(kmerRegex.FindStringSubmatch(field)[1])
			prof := NewKmerProfile(k)
			prof.GetSeqKmers(seq)
			prof.CountsToFreqs()
			return prof.KmerDist(ref[k]), err
		}
		err = fmt.Errorf("Invalid metric: %s", field)
	}
	return math.NaN(), err
}
