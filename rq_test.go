package qr

import (
	"math"
	"testing"
)

func TestRQ(t *testing.T) {
	A := [][]float64{
		[]float64{1, 3, -1},
		[]float64{2, -0.5, 2.2},
		[]float64{3, -0.4, 2.1},
		[]float64{3, 2.1, 1.5},
	}
	x := []float64{1.2, 2.1, 3.4}
	b := matvec(A, x)

	d, err := NewReal(A)
	check(t, err)

	r, err := d.Solve(b)
	check(t, err)
	if len(r) != len(x) {
		t.Fatal()
	}
	for i := range r {
		if e := math.Abs(r[i] - x[i]); e > 3e-15 {
			t.Fatalf("expected %v got %v e=%v\n", x[i], r[i], e)
		}
	}
}
func check(t *testing.T, e error) {
	if e != nil {
		t.Fatal(e)
	}
}

func matvec(A [][]float64, v []float64) []float64 {
	r := make([]float64, len(A))
	for i := range r {
		s := 0.0
		for k := range v {
			s += A[i][k] * v[k]
		}
		r[i] = s
	}
	return r
}
