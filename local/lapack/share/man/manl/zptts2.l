.TH  ZPTTS2 1 "April 2011" " LAPACK routine (version 3.3.1) " " LAPACK routine (version 3.3.1) " 
.SH NAME
LAPACK-3 - solves a tridiagonal system of the form  A * X = B using the factorization A = U**H *D*U or A = L*D*L**H computed by ZPTTRF
.SH SYNOPSIS
.TP 19
SUBROUTINE ZPTTS2(
IUPLO, N, NRHS, D, E, B, LDB )
.TP 19
.ti +4
INTEGER
IUPLO, LDB, N, NRHS
.TP 19
.ti +4
DOUBLE
PRECISION D( * )
.TP 19
.ti +4
COMPLEX*16
B( LDB, * ), E( * )
.SH PURPOSE
ZPTTS2 solves a tridiagonal system of the form
   A * X = B
using the factorization A = U**H *D*U or A = L*D*L**H computed by ZPTTRF.
 D is a diagonal matrix specified in the vector D, U (or L) is a unit
 bidiagonal matrix whose superdiagonal (subdiagonal) is specified in
 the vector E, and X and B are N by NRHS matrices.
.br
.SH ARGUMENTS
.TP 8
 IUPLO   (input) INTEGER
 Specifies the form of the factorization and whether the
 vector E is the superdiagonal of the upper bidiagonal factor
 U or the subdiagonal of the lower bidiagonal factor L.
 = 1:  A = U**H *D*U, E is the superdiagonal of U
.br
 = 0:  A = L*D*L**H, E is the subdiagonal of L
.TP 8
 N       (input) INTEGER
 The order of the tridiagonal matrix A.  N >= 0.
.TP 8
 NRHS    (input) INTEGER
 The number of right hand sides, i.e., the number of columns
 of the matrix B.  NRHS >= 0.
.TP 8
 D       (input) DOUBLE PRECISION array, dimension (N)
 The n diagonal elements of the diagonal matrix D from the
 factorization A = U**H *D*U or A = L*D*L**H.
.TP 8
 E       (input) COMPLEX*16 array, dimension (N-1)
 If IUPLO = 1, the (n-1) superdiagonal elements of the unit
 bidiagonal factor U from the factorization A = U**H*D*U.
 If IUPLO = 0, the (n-1) subdiagonal elements of the unit
 bidiagonal factor L from the factorization A = L*D*L**H.
.TP 8
 B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
 On entry, the right hand side vectors B for the system of
 linear equations.
 On exit, the solution vectors, X.
.TP 8
 LDB     (input) INTEGER
 The leading dimension of the array B.  LDB >= max(1,N).