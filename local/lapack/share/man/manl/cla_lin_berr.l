.TH  CLA_LIN_BERR 1 "April 2011" "    LAPACK routine (version 3.2.2)                                 " "    LAPACK routine (version 3.2.2)                                 " 
.SH NAME
LAPACK-3 - CLA_LIN_BERR compute componentwise relative backward error from  the formula  max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )  where abs(Z) is the componentwise absolute value of the matrix  or vector Z
.SH SYNOPSIS
.TP 24
SUBROUTINE CLA_LIN_BERR
( N, NZ, NRHS, RES, AYB, BERR )
.TP 24
.ti +4
IMPLICIT
NONE
.TP 24
.ti +4
INTEGER
N, NZ, NRHS
.TP 24
.ti +4
REAL
AYB( N, NRHS ), BERR( NRHS )
.TP 24
.ti +4
COMPLEX
RES( N, NRHS )
.SH PURPOSE
   CLA_LIN_BERR computes componentwise relative backward error from
   the formula
       max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
   where abs(Z) is the componentwise absolute value of the matrix
   or vector Z.
    N       (input) INTEGER
.br
    The number of linear equations, i.e., the order of the
    matrix A.  N >= 0.
.br
    NZ      (input) INTEGER
.br
    We add (NZ+1)*SLAMCH( \(aqSafe minimum\(aq ) to R(i) in the numerator to
    guard against spuriously zero residuals. Default value is N.
    NRHS    (input) INTEGER
.br
    The number of right hand sides, i.e., the number of columns
    of the matrices AYB, RES, and BERR.  NRHS >= 0.
.br
    RES    (input) DOUBLE PRECISION array, dimension (N,NRHS)
    The residual matrix, i.e., the matrix R in the relative backward
    error formula above.
.br
    AYB    (input) DOUBLE PRECISION array, dimension (N, NRHS)
    The denominator in the relative backward error formula above, i.e.,
    the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B
    are from iterative refinement (see cla_gerfsx_extended.f).
    
.br
    BERR   (output) COMPLEX array, dimension (NRHS)
.br
    The componentwise relative backward error from the formula above.