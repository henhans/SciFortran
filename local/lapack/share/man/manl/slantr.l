.TH  SLANTR 1 "April 2011" " LAPACK auxiliary routine (version 3.2) " " LAPACK auxiliary routine (version 3.2) " 
.SH NAME
LAPACK-3 - returns the value of the one norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a trapezoidal or triangular matrix A
.SH SYNOPSIS
.TP 14
REAL FUNCTION
SLANTR( NORM, UPLO, DIAG, M, N, A, LDA,
WORK )
.TP 14
.ti +4
CHARACTER
DIAG, NORM, UPLO
.TP 14
.ti +4
INTEGER
LDA, M, N
.TP 14
.ti +4
REAL
A( LDA, * ), WORK( * )
.SH PURPOSE
SLANTR  returns the value of the one norm,  or the Frobenius norm, or
the  infinity norm,  or the  element of  largest absolute value  of a
trapezoidal or triangular matrix A.
.SH DESCRIPTION
 SLANTR returns the value
.br
    SLANTR = ( max(abs(A(i,j))), NORM = \(aqM\(aq or \(aqm\(aq
.br
             (
.br
             ( norm1(A),         NORM = \(aq1\(aq, \(aqO\(aq or \(aqo\(aq
.br
             (
.br
             ( normI(A),         NORM = \(aqI\(aq or \(aqi\(aq
.br
             (
.br
             ( normF(A),         NORM = \(aqF\(aq, \(aqf\(aq, \(aqE\(aq or \(aqe\(aq
 where  norm1  denotes the  one norm of a matrix (maximum column sum),
 normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
 normF  denotes the  Frobenius norm of a matrix (square root of sum of
 squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
.SH ARGUMENTS
.TP 8
 NORM    (input) CHARACTER*1
 Specifies the value to be returned in SLANTR as described
 above.
.TP 8
 UPLO    (input) CHARACTER*1
 Specifies whether the matrix A is upper or lower trapezoidal.
 = \(aqU\(aq:  Upper trapezoidal
.br
 = \(aqL\(aq:  Lower trapezoidal
 Note that A is triangular instead of trapezoidal if M = N.
.TP 8
 DIAG    (input) CHARACTER*1
 Specifies whether or not the matrix A has unit diagonal.
 = \(aqN\(aq:  Non-unit diagonal
.br
 = \(aqU\(aq:  Unit diagonal
.TP 8
 M       (input) INTEGER
 The number of rows of the matrix A.  M >= 0, and if
 UPLO = \(aqU\(aq, M <= N.  When M = 0, SLANTR is set to zero.
.TP 8
 N       (input) INTEGER
 The number of columns of the matrix A.  N >= 0, and if
 UPLO = \(aqL\(aq, N <= M.  When N = 0, SLANTR is set to zero.
.TP 8
 A       (input) REAL array, dimension (LDA,N)
 The trapezoidal matrix A (A is triangular if M = N).
 If UPLO = \(aqU\(aq, the leading m by n upper trapezoidal part of
 the array A contains the upper trapezoidal matrix, and the
 strictly lower triangular part of A is not referenced.
 If UPLO = \(aqL\(aq, the leading m by n lower trapezoidal part of
 the array A contains the lower trapezoidal matrix, and the
 strictly upper triangular part of A is not referenced.  Note
 that when DIAG = \(aqU\(aq, the diagonal elements of A are not
 referenced and are assumed to be one.
.TP 8
 LDA     (input) INTEGER
 The leading dimension of the array A.  LDA >= max(M,1).
.TP 8
 WORK    (workspace) REAL array, dimension (MAX(1,LWORK)),
 where LWORK >= M when NORM = \(aqI\(aq; otherwise, WORK is not
 referenced.