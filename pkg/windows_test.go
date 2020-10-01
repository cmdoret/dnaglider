package pkg

import (
	"fmt"
	"testing"
)

func TestChunkGenome(t *testing.T) {
	var pos, nChunks int
	var prevID string
	records := StreamGenome("../tests/genome.fa", 1)
	chunks := ChunkGenome(records, 5, 10)
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
