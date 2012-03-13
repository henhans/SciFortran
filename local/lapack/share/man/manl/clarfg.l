.TH  CLARFG 1 "April 2011" " LAPACK auxiliary routine (version 3.3.1) " " LAPACK auxiliary routine (version 3.3.1) " 
.SH NAME
LAPACK-3 - generates a complex elementary reflector H of order n, such that   H**H * ( alpha ) = ( beta ), H**H * H = I
.SH SYNOPSIS
.TP 19
SUBROUTINE CLARFG(
N, ALPHA, X, INCX, TAU )
.TP 19
.ti +4
INTEGER
INCX, N
.TP 19
.ti +4
COMPLEX
ALPHA, TAU
.TP 19
.ti +4
COMPLEX
X( * )
.SH PURPOSE
CLARFG generates a complex elementary reflector H of order n, such
that
              (   x   )   (   0  )
.br
 where alpha and beta are scalars, with beta real, and x is an
 (n-1)-element complex vector. H is represented in the form
       H = I - tau * ( 1 ) * ( 1 v**H ) ,
.br
                     ( v )
.br
 where tau is a complex scalar and v is a complex (n-1)-element
 vector. Note that H is not hermitian.
.br
 If the elements of x are all zero and alpha is real, then tau = 0
 and H is taken to be the unit matrix.
.br
 Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
.br
.SH ARGUMENTS
.TP 8
 N       (input) INTEGER
 The order of the elementary reflector.
.TP 8
 ALPHA   (input/output) COMPLEX
 On entry, the value alpha.
 On exit, it is overwritten with the value beta.
.TP 8
 X       (input/output) COMPLEX array, dimension
 (1+(N-2)*abs(INCX))
 On entry, the vector x.
 On exit, it is overwritten with the vector v.
.TP 8
 INCX    (input) INTEGER
 The increment between elements of X. INCX > 0.
.TP 8
 TAU     (output) COMPLEX
 The value tau.