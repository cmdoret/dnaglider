package pkg

import (
	"fmt"
	"strconv"
	"testing"
)

var (
	CHUNKSIZE int = 10
	WINSIZE   int = 5
)

func TestChunkGenome(t *testing.T) {
	var pos, nChunks int
	var prevID string
	records := StreamGenome("../tests/genome.fa", 1)
	chunks := ChunkGenome(records, WINSIZE, CHUNKSIZE)
	for chunk := range chunks {
		nChunks++
		if len(chunk.Seq.Seq) != 50 {
			fmt.Println(len(chunk.Seq.Seq))
			t.Errorf("Sequence length of chunk is not wsize * chunksize")
		}
		if string(chunk.ID) == prevID {
			if pos != chunk.BpStart {
				t.Errorf("Chunk does not start at the end of previous chunk")
			}
			pos = chunk.BpEnd
		} else {
			pos = 0
		}
	}

}
func TestConsumeChunks(t *testing.T) {
	records := StreamGenome("../tests/genome.fa", 1)
	ref := make(map[int]KmerProfile)
	ref[3] = FastaToKmers("../tests/genome.fa", 3)
	chunks := ChunkGenome(records, WINSIZE, CHUNKSIZE)
	expHeader := []string{"chrom", "start", "end", "GC", "GCSKEW", "ENTRO"}
	results := ConsumeChunks(chunks, []string{"GC", "GCSKEW", "ENTRO"}, ref)
	for res := range results {
		for i := range res.Header {
			if res.Header[i] != expHeader[i] {
				t.Errorf("Result header is incorrect.: %s", res.Header)
			}
		}
		if len(res.Data) != CHUNKSIZE {
			t.Errorf(
				"Chunk result contains incorrect number of windows: %d instead of %d",
				len(res.Data),
				CHUNKSIZE,
			)
		}
		winStart, _ := strconv.Atoi(res.Data[0][1])
		winEnd, _ := strconv.Atoi(res.Data[0][2])
		if (winEnd-winStart)+1 != WINSIZE {
			t.Errorf(
				"First window has incorrect reported length in chunk result: %d instead of %d",
				(winEnd-winStart)+1,
				WINSIZE,
			)
		}
	}

}

func TestMakeRange(t *testing.T) {
	var start, end, step int
	start = 3
	end = 12
	step = 3
	var tRange []int
	tRange = MakeRange(start, end, step)
	obsLen := len(tRange)
	expLen := 1 + (end-start)/step
	if obsLen != expLen {
		t.Errorf("Length of MakeRange output is incorrect (%d instead of %d)", obsLen, expLen)
	}
	if tRange[obsLen-1] != end {
		t.Errorf("Last element of MakeRange output is incorrect")
	}
	if tRange[0] != start {
		t.Errorf("First element of MakeRange output is incorrect")
	}
}

func TestMinInt(t *testing.T) {
	var a, b, m int
	a = 3
	b = 4
	m = MinInt(a, b)
	if m != a {
		t.Errorf("MinInt returned wrong min value")
	}
}

func TestBuild2dSlice(t *testing.T) {
	rows := 11
	cols := 8
	slice2d := Build2dSlice(rows, cols)
	if len(slice2d) != rows {
		t.Errorf("Wrong number of rows: %d instead of %d", len(slice2d), rows)
	}
	if len(slice2d[0]) != cols {
		t.Errorf("Wrong number of cols: %d instead of %d", len(slice2d[0]), cols)
	}
	for idx, row := range slice2d {
		if len(row) != cols {
			t.Errorf("2D slice is not rectangle: wrong number of elements at row %d", idx)
		}
	}
}
