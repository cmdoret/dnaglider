package pkg

import (
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

// ChunkResult stores the computed statistics of windows from a chunk in the
// form of a table. Each row is a window, each column is a feature.
type ChunkResult struct {
	Header []string
	Data   [][]string
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

// Build2dSlice builds a 2d slice of float64 of target size
func Build2dSlice(rows int, cols int) [][]string {
	slice2d := make([][]string, rows)
	for i := range slice2d {
		slice2d[i] = make([]string, cols)
	}
	return slice2d
}

// ChunkGenome receives fastx.Record from a channel and produces a channel of record
// chunks. Each chunk contains multiple windows. The chunkSize is given in number of
// windows, and the windows size is in basepair
func ChunkGenome(records <-chan fastx.Record, winSize int, chunkSize int) <-chan Chunk {
	chunkLen := chunkSize * winSze
	chunks := make(chan Chunk, 3)
	var bpStart, bpEnd int
	var bpStart, bpEnd int
	go func() {
		for rec := range records {
			seqLen = len(rec.Seq.Seq)
			bpEnd = 0
			bpStart = 1
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

// ConsumeChunks computes window-based statistics in chunks and stores them in a ChunkResult struct.
func ConsumeChunks(chunks <-chan Chunk, metrics []string, refProfile [int]KmerProfile) chan ChunkResult {
	// There are 3 columns for coordinates (chrom start end), and 1 per feature
	nFeatures := 3 + len(refProfile) + len(metrics)
	// Generate column names
	header := []string{"chrom", "start", "end"}
	header = append(header, metrics)
	header = append()
	for i, col := range  {
		header[i] = col
	}

	go func() {
		for chunk := range chunks {
			nWindows := len(chunk.Starts)
			results := ChunkResult{header, Build2dSlice(nWindows, nFeatures)}
			for i, start := range chunk.Starts {
				result.Data[i][0] = chunk.ID
				result.Data[i][1] = chunk.BpStart + start
				result.Data[i][2] = chunk.BpStart + start + chunk.wSize
			}
		}
		close(out)
	}()
	return out
}
