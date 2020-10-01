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
	var bpStart, bpEnd int
	chunkLen := chunkSize * winSize
	wins := make(chan Chunk)
	go func() {
		// Need to add overlap between chunks
		for rec := range records {
			seqLen := len(rec.Seq.Seq)
			bpEnd = MinInt(bpStart+chunkLen, seqLen)
			chunk := Chunk{
				ID:      rec.ID,
				BpStart: bpStart,
				BpEnd:   bpEnd,
				Starts:  MakeRange(0, seqLen, winSize),
				wSize:   winSize,
			}
			wins <- chunk
			bpStart = bpEnd
		}
		close(wins)
	}()
	return wins
}

// ConsumeChunks computes window-based statistics in chunks
func ConsumeChunks(chunks <-chan Chunk) {
	var end int
	go func() {
		for chunk := range chunks {
			winStart := chunk.BpStart
			for start := range chunk.Starts {
				end = start + chunk.wSize
				gc := SeqGC(chunk.Seq.SubSeq(start, end))
				fmt.Printf("%s\t%d\t%d\t%f", chunk.ID, winStart, winStart+chunk.wSize, gc)
				winStart += chunk.wSize

			}
		}
	}()
}
