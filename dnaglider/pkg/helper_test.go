package pkg

import (
	"testing"
)

func TestCheckFields(t *testing.T) {
	validFields := []string{"GC", "GCSKEW", "ATSKEW", "ENTRO", "3MER", "14MER"}
	invalidFields := []string{"John Smith", "paper"}
	if err := checkFields(validFields); err != nil {
		t.Errorf("Valid fields raised an error")
	}
	if err := checkFields(invalidFields); err == nil {
		t.Errorf("Invalid fields did not raise an error")
	}
}
