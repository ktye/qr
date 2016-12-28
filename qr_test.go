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

func ExampleSolve() {
	// Solution to the overdetermined 3x2 system: A x = b.
	// Matrix A (3x2).
	A := [][]complex128{
		[]complex128{complex(1, 2), complex(5, -1)},
		[]complex128{complex(6, 8), complex(1, 0)},
		[]complex128{complex(3, -2), complex(-7, 3)},
	}

	// Right hand side (3x1):
	b := []complex128{
		complex(5, -2),
		complex(-3, 1),
		complex(3, 0),
	}

	// Compute the QR decomposition.
	if D, err := New(A); err != nil {
		panic(err)
	} else {
		// Solve the system.
		if x, err := D.Solve(b); err != nil {
			panic(err)
		} else {
			fmt.Println(x)
		}
	}
	// Output: [(0.040087623220153294+0.20569550930996708i) (0.10186199342825833-0.1207009857612268i)]
}
