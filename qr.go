// Package qr does a QR decomposition of a general complex128 matrix.
//
// Reference:
// Walter Gander, Martin J. Gander, Felix Kwok:
// Scientific Computing, an Introduction Using Maple and Matlab,
// Springer, April 2014, ISBN 978-3-319-04324-1
// Pages 359-361
package qr

import (
	"errors"
	"math"
	"math/cmplx"
)

// QR is the result of a QR decomposition in a compact storage.
// Householder vectors are stored in the lower part of the matrix.
// H[k][k:] is the kth Householder vector.
// In the upper part the matrix R is stored without the diagonal.
type QR struct {
	H     [][]complex128
	Rdiag []complex128 // Missing diagonal of R.
	m, n  int          // number of rows and columns
}

// New calculates the QR Decomposition of a rectangular matrix.
func New(A [][]complex128) (QR, error) {
	m := len(A)    // Number of rows.
	n := len(A[0]) // Number of columns.
	if m < n {
		return QR{}, errors.New("qr: matrix is underdetermined")
	}

	// Build the workspace with a working copy of A in colum major format.
	H := make([][]complex128, n)
	Rdiag := make([]complex128, n)
	for i := 0; i < n; i++ {
		H[i] = make([]complex128, m)
		for k := 0; k < m; k++ {
			H[i][k] = A[k][i]
		}
	}

	for j := 0; j < n; j++ {
		s := VectorNorm(H[j][j:])
		if s == 0 {
			return QR{}, errors.New("matrix contains zero-columns")
		}

		Rdiag[j] = -complex(s, 0) * cmplx.Rect(1, cmplx.Phase(H[j][j])) // Diagonal element.
		f := complex(math.Sqrt(s*(s+cmplx.Abs(H[j][j]))), 0)
		H[j][j] -= Rdiag[j]

		for k := j; k < m; k++ {
			H[j][k] /= f
		}
		for i := j + 1; i < n; i++ {
			var sum complex128
			for k := j; k < m; k++ {
				sum += cmplx.Conj(H[j][k]) * H[i][k]
			}
			for k := j; k < m; k++ {
				H[i][k] -= H[j][k] * sum
			}
		}
	}

	return QR{
		H:     H,
		Rdiag: Rdiag,
		m:     m,
		n:     n,
	}, nil
}

// Solve solves the overdetermined system A*x = b.
func (D QR) Solve(b []complex128) ([]complex128, error) {
	if len(b) != D.m {
		return nil, errors.New("qr: wrong input dimension for QR.Solve.")
	}
	if QTx, err := D.QMul(b); err != nil {
		return nil, err
	} else {
		return D.RSolve(QTx)
	}
}

// LeastSquareSolve solves the overdetermined system A*x = b.
func LeastSquareSolve(A [][]complex128, b []complex128) ([]complex128, error) {
	if D, err := New(A); err != nil {
		return nil, err
	} else {
		return D.Solve(b)
	}
}

// QMul does a matrix vector multiplication of matrix Q' from a QR decomposition with vector x.
func (D QR) QMul(x []complex128) ([]complex128, error) {
	if len(x) != D.m {
		return nil, errors.New("qr: input vector lengths mismatch for QMul.")
	}
	y := make([]complex128, D.m)
	for i := 0; i < D.m; i++ {
		y[i] = x[i]
	}
	for j := 0; j < D.n; j++ {
		var sum complex128
		for k := j; k < D.m; k++ {
			sum += cmplx.Conj(D.H[j][k]) * y[k]
		}
		for k := j; k < D.m; k++ {
			y[k] -= D.H[j][k] * sum
		}
	}
	return y, nil
}

// RSolve solves the system R*x = b with R of the QR decomposition using back-substitution.
func (D QR) RSolve(b []complex128) ([]complex128, error) {
	if len(b) != D.m {
		return nil, errors.New("qr: input vector lengths mismatch for RSolve.")
	}
	x := make([]complex128, D.m)
	for i := 0; i < D.m; i++ {
		x[i] = b[i]
	}
	for i := D.n - 1; i >= 0; i-- {
		x[i] = b[i]
		for j := i + 1; j < D.n; j++ {
			x[i] -= D.H[j][i] * x[j]
		}
		x[i] /= D.Rdiag[i]
	}
	return x[0:D.n], nil
}

// VectorNorm computes the vector norm of a complex128 vector without over/underflow.
func VectorNorm(x []complex128) (norm float64) {
	for i := 0; i < len(x); i++ {
		norm = math.Hypot(norm, cmplx.Abs(x[i]))
	}
	return
}
