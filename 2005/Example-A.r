#  
# Exponential Distribution Model
#
# Ref 1 = Usher and Hodgson (1988), IEEE Reliability Vol.37 pg. 550 - 555 
#
#=========================================================
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2005/Rs3.R")
#=========================================================
X1 = c(0.021, 0.038, 0.054, 0.066, 0.076, 0.078, 0.123, 0.130, 0.152, 0.159,
       0.199, 0.201, 0.204, 0.215, 0.218, 0.281, 0.295, 0.310, 0.338, 0.341,
       0.354, 0.358, 0.431, 0.457, 0.545, 0.569, 0.677, 0.818, 0.946, 1.486)
#===============================================================================
maxit = 10; outputdigit=3; 

## Complete Information (No Masking)
cat("\n=======================\n")
cat(" Complete Information")
cat("\n-----------------------\n")
modes = list(2, 2, 3, 3, 2, 2, 3, 3, 2, 1,
       3, 1, 1, 3, 2, 1, 2, 3, 3, 2,
       1, 2, 1, 3, 1, 2, 3, 2, 2, 1)
OUTP = expo.cm.EM1(X1, modes, start=1, maxits=maxit, eps=0.0001)
OUTAB= expo.cm.EM2(X1, modes, start=1, maxits=maxit, eps=0.0001)
OUTPUT = cbind( OUTAB$para, OUTP$para )
colnames(OUTPUT) = c("AB", "P")
rownames(OUTPUT) = c("lambda1", "lambda2", "lambda3")
round( OUTPUT, outputdigit )

## Partial Masking  (General Case. Page 554 of Ref 1)
cat("\n=======================\n")
cat(" General Masking")
cat("\n-----------------------\n")
modes = list(2,1:2,  3,  3,1:2,2:3, 3,c(1,3),1:3, 1,
          3,  1,  1,2:3,1:2,  1, 2,     3,  3, 2,
          1,  2,1:3,  3,1:3,  2, 3,     2,2:3, 1)
OUTP = expo.cm.EM1(X1, modes, start=1, maxits=maxit, eps=0.0001)
OUTAB= expo.cm.EM2(X1, modes, start=1, maxits=maxit, eps=0.0001)
OUTPUT = cbind( OUTAB$para, OUTP$para )
colnames(OUTPUT) = c("AB", "P")
rownames(OUTPUT) = c("lambda1", "lambda2", "lambda3")
round( OUTPUT, outputdigit )

## Compelte Masking  (Case 1. Page 554 of Ref 1)
cat("\n=======================\n")
cat(" Case 1 (Complete Masking)")
cat("\n-----------------------\n")
modes = list(2, 2, 3,  3, 2, 2, 3, 3, 1:3, 1,
          3, 1, 1,  3, 2, 1, 2, 3, 3, 2,
          1, 2,1:3, 3,1:3,2, 3, 2, 2, 1)
OUTP = expo.cm.EM1(X1, modes, start=1, maxits=maxit)
OUTAB= expo.cm.EM2(X1, modes, start=1, maxits=maxit)
OUTPUT = cbind( OUTAB$para, OUTP$para )
colnames(OUTPUT) = c("AB", "P")
rownames(OUTPUT) = c("lambda1", "lambda2", "lambda3")
round( OUTPUT, outputdigit )


## Partial Masking  (Case 2. Page 554 of Ref 1)
cat("\n=======================\n")
cat(" Case 2")
cat("\n-----------------------\n")
modes = list(2,1:2, 3,  3,1:2, 2, 3, 3, 2, 1,
          3,  1, 1,  3,1:2, 1, 2, 3, 3, 2,
          1,  2, 1,  3,  1, 2, 3, 2, 2, 1)
OUTP = expo.cm.EM1(X1, modes, start=1, maxits=maxit)
OUTAB= expo.cm.EM2(X1, modes, start=1, maxits=maxit)
OUTPUT = cbind( OUTAB$para, OUTP$para )
colnames(OUTPUT) = c("AB", "P")
rownames(OUTPUT) = c("lambda1", "lambda2", "lambda3")
round( OUTPUT, outputdigit )

## Partial Masking  (Case 3. Page 554 of Ref 1)
cat("\n=======================\n")
cat(" Case 3")
cat("\n-----------------------\n")
modes = list(2,1:2,  3,  3,1:2, 2, 3, 3,1:3, 1,
          3,  1,  1,  3,1:2, 1, 2, 3,  3, 2,
          1,  2,1:3,  3,1:3, 2, 3, 2,  2, 1)
OUTP = expo.cm.EM1(X1, modes, start=1, maxits=maxit)
OUTAB= expo.cm.EM2(X1, modes, start=1, maxits=maxit)
OUTPUT = cbind( OUTAB$para, OUTP$para )
colnames(OUTPUT) = c("AB", "P")
rownames(OUTPUT) = c("lambda1", "lambda2", "lambda3")
round( OUTPUT, outputdigit )

