# qr
QR decomposition of a complex matrix and least squares solver

[![GoDoc](https://godoc.org/github.com/ktye/qr?status.svg)](https://godoc.org/github.com/ktye/qr)

The package provides an implementation of the QR decomposition for general complex matrices for the go programming language.
It uses the Householder algorithm described in this book:

> Walter Gander, Martin J. Gander, Felix Kwok:
> Scientific Computing, an Introduction Using Maple and Matlab,
> Springer, April 2014, ISBN 978-3-319-04324-1
> Pages 359-361.

The QR decomposition is often used to solve overdetermined linear systems of equations in the form
> A x = b

where A is a complex matrix with size m x n and m >= n.
b is the right hand side vector of length m and x the vector of unknowns with the same size as b.

The package is self-contained and uses only the standard library.

Example:
```go
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
                        // do something with x
                }
        }

```
