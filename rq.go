package qr

import (
	"errors"
	"math"
)

type RQ struct {
	H     [][]float64
	Rdiag []float64
	m, n  int
}

func NewReal(A [][]float64) (RQ, error) {
	m := len(A)    // Number of rows.
	n := len(A[0]) // Number of columns.
	if m < n {
		return RQ{}, errors.New("qr: matrix is underdetermined")
	}
	H := make([][]float64, n)
	Rdiag := make([]float64, n)
	for i := 0; i < n; i++ {
		H[i] = make([]float64, m)
		for k := 0; k < m; k++ {
			H[i][k] = A[k][i]
		}
	}
	for j := 0; j < n; j++ {
		s := RealNorm(H[j][j:])
		if s == 0 {
			return RQ{}, errors.New("matrix contains zero-columns")
		}
		if H[j][j] > 0 {
			Rdiag[j] = -s
		} else {
			Rdiag[j] = s
		}
		f := 1.0 / math.Sqrt(s*(s+math.Abs(H[j][j])))
		H[j][j] -= Rdiag[j]
		for k := j; k < m; k++ {
			H[j][k] *= f
		}
		for i := j + 1; i < n; i++ {
			var sum float64
			for k := j; k < m; k++ {
				sum += H[j][k] * H[i][k]
			}
			for k := j; k < m; k++ {
				H[i][k] -= H[j][k] * sum
			}
		}
	}
	return RQ{
		H:     H,
		Rdiag: Rdiag,
		m:     m,
		n:     n,
	}, nil
}
func (D RQ) Solve(b []float64) ([]float64, error) {
	if len(b) != D.m {
		return nil, errors.New("qr: wrong input dimension for QR.Solve.")
	}
	if QTx, err := D.QMul(b); err != nil {
		return nil, err
	} else {
		return D.RSolve(QTx)
	}
}
func (D RQ) QMul(x []float64) ([]float64, error) {
	if len(x) != D.m {
		return nil, errors.New("qr: input vector lengths mismatch for QMul.")
	}
	y := make([]float64, D.m)
	for i := 0; i < D.m; i++ {
		y[i] = x[i]
	}
	for j := 0; j < D.n; j++ {
		var sum float64
		for k := j; k < D.m; k++ {
			sum += D.H[j][k] * y[k]
		}
		for k := j; k < D.m; k++ {
			y[k] -= D.H[j][k] * sum
		}
	}
	return y, nil
}
func (D RQ) RSolve(b []float64) ([]float64, error) {
	if len(b) != D.m {
		return nil, errors.New("qr: input vector lengths mismatch for RSolve.")
	}
	x := make([]float64, D.m)
	for i := 0; i < D.m; i++ {
		x[i] = b[i]
	}
	for i := D.n - 1; i >= 0; i-- {
		s := 0.0
		for j := i + 1; j < D.n; j++ {
			s += D.H[j][i] * x[j]
			x[i] -= D.H[j][i] * x[j]
		}
		x[i] /= D.Rdiag[i]
	}
	return x[0:D.n], nil
}
func RealNorm(v []float64) (r float64) {
	s := 0.0
	for _, x := range v {
		if x != 0 {
			x = math.Abs(x)
			if s < x {
				t := s / x
				r = 1 + r*t*t
				s = x
			} else {
				t := x / s
				r += t * t
			}
		}
	}
	return s * math.Sqrt(r)
}
