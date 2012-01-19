cc -O sine_gen.c -lm -o bin/sine_gen
cc -O rand_gen.c -lm -o bin/rand_gen
cc -O scale.c 	 -lm -o bin/scale
cc -O shift.c 	 -lm -o bin/shift
cc -O normalize.c -lm -o bin/normalize
cc -O cmag.c 	 -lm -o bin/cmag
cc -O db.c 	 -lm -o bin/db
cc -O phase.c 	 -lm -o bin/phase
cc -O real2cmplx.c -lm -o bin/real2cmplx
cc -O extract_real.c	  -lm -o bin/extract_real
cc -O extract_imaginary.c -lm -o bin/extract_imaginary
cc -O add.c 	 -lm -o bin/add
cc -O sub.c 	 -lm -o bin/sub
cc -O conjg.c 	 -lm -o bin/conjg
cc -O dotpr.c 	 -lm -o bin/dotpr
cc -O cdotpr.c 	 -lm -o bin/cdotpr
cc -O matmul.c	 -lm -o bin/matmul
cc -O multpbp.c	 -lm -o bin/multpbp
cc -O outpr.c	 -lm -o bin/outpr
cc -O zfill.c	 -lm -o bin/zfill
cc -O statistics.c -lm -o bin/statistics
cc -O ucirc.c	 -lm -o bin/ucirc
cc -O window.c	 -lm -o bin/window
cc -O transpose.c -lm -o bin/transpose
cc -O fft.c	 -lm -o bin/fft
cc -O convolve.c -lm -o bin/convolve
cc -O correlate.c -lm -o bin/correlate
cc -O fft2D.c	 -lm -o bin/fft2D
cc -O matXvec.c	 -lm -o bin/matXvec
cc -O matinvert.c -lm -o bin/matinvert
cc -O mdeterm.c	 -lm -o bin/mdeterm
cc -O augmatvec.c -lm -o bin/augmatvec
cc -O join2xy.c	 -lm -o bin/join2xy
cc -O mat2image.c -lm -o bin/mat2image
cc -O histogram.c -lm -o bin/histogram
cc -O barchart.c -lm -o bin/barchart
cc -O differentiate.c -lm -o bin/differentiate
cc -O integrate.c -lm -o bin/integrate
cc -O nums2vec.c -lm -o bin/nums2vec
cc -O vec2nums.c -lm -o bin/vec2nums
cc -O enter_vec.c -lm -o bin/enter_vec
cc -O abs.c -lm -o bin/abs
cc -O image2matrix.c -lm -o bin/image2matrix
cc -O decimate.c -lm -o bin/decimate
