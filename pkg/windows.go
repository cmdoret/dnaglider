package pkg

import (
	"fmt"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

// Chunk is a piece of fastx.Record sequence, containing the associated
// id genomic coordinate. It also contains indices of sliding windows
// in which statistics will be computed.
type Chunk struct {
	ID      []byte
	BpStart int
	BpEnd   int
	wSize   int
	Seq     *seq.Seq
	Starts  []int
}

// MinInt returns the smallest of two integers. If both are equal, the second
// input is returned.
func MinInt(x, y int) int {
	if x < y {
		return x
	}
	return y
}

// MakeRange generates a slice of ints from start to end, where each value is
// spaced by step
func MakeRange(start, end, step int) []int {
	arr := make([]int, 0, 1+(end-start)/step)
	for start <= end {
		arr = append(arr, start)
		start += step
	}
	return arr
}

// ChunkGenome receives fastx.Record from a channel and produces a channel of record
// chunks. Each chunk contains multiple windows. The chunkSize is given in number of
// windows, and the windows size is in basepair
func ChunkGenome(records <-chan fastx.Record, winSize int, chunkSize int) <-chan Chunk {
	chunkLen := chunkSize * winSize
	chunks := make(chan Chunk, 3)
	go func() {
		var bpStart, bpEnd int
		// Need to add overlap between chunks
		for rec := range records {
			seqLen := len(rec.Seq.Seq)
			bpStart = 1
			bpEnd = 0
			for bpEnd < seqLen {
				bpEnd = MinInt(bpStart+chunkLen, seqLen)
				chunk := Chunk{
					ID:      rec.ID,
					BpStart: bpStart,
					BpEnd:   bpEnd,
					Starts:  MakeRange(0, bpEnd-bpStart, winSize),
					wSize:   winSize,
					Seq:     rec.Seq.SubSeq(bpStart, bpEnd),
				}

				chunks <- chunk
				bpStart = bpEnd
			}
		}
		close(chunks)
		return
	}()
	return chunks
}

// ConsumeChunks computes window-based statistics in chunks
func ConsumeChunks(chunks <-chan Chunk) chan string {
	var chunkStr string
	out := make(chan string, 1)
	go func() {
		for chunk := range chunks {
			chunkStr = ""
			for _, start := range chunk.Starts {
				chunkStr += fmt.Sprintf(
					"%s\t%d\t%d\t%f\n",
					chunk.ID,
					chunk.BpStart+start,
					chunk.BpStart+start+chunk.wSize,
					SeqGC(chunk.Seq.SubSeq(start+1, start+chunk.wSize)),
				)
			}
			out <- chunkStr
		}
		close(out)
		return
	}()
	return out
}
