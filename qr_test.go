package qr

import (
	"fmt"
	"math/cmplx"
	"testing"
)

func checkVectorSimilarity(a, b []complex128) error {
	if len(a) != len(b) {
		return fmt.Errorf("dimensions mismatch %d != %d", len(a), len(b))
	}
	eps := 1E-14
	for i := 0; i < len(a); i++ {
		if diff := cmplx.Abs(a[i] - b[i]); diff > eps {
			return fmt.Errorf("element %d differs by %v", i, diff)
		}
	}
	return nil
}

func TestLeastSquareSolve(t *testing.T) {
	A := lsqrA
	b := lsqrB
	X := lsqrX
	D, err := New(A)
	if err != nil {
		t.Error()
	}
	if x, err := D.Solve(b); err != nil {
		t.Fatal(err)
	} else {
		if err := checkVectorSimilarity(x, X); err != nil {
			t.Error(err)
		}
	}
}
