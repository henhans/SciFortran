.TH  SSCAL 1 "April 2011" "BLAS routine" "BLAS routine" 
.SH NAME
LAPACK-3 - scales a vector by a constant
.SH SYNOPSIS
.TP 31
SUBROUTINE SSCAL(N,SA,SX,INCX)
.TP 31
.ti +4
REAL
SA
.TP 31
.ti +4
INTEGER
INCX,N
.TP 31
.ti +4
REAL
SX(*)
.SH PURPOSE
   scales a vector by a constant.
    uses unrolled loops for increment equal to 1.
.br

.SH FURTHER DETAILS
    jack dongarra, linpack, 3/11/78.
.br
    modified 3/93 to return if incx .le. 0.
.br
    modified 12/3/93, array(1) declarations changed to array(*)
