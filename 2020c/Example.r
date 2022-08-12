# Testing H0: C <= 1.33 vs. H1: C > 1.33

conf = 0.8 # <- Type confidence level

# It needs boot R package. 
# If boot package is not installed, type
# install.packages("boot")
library(boot)


############################################################
data = c(
5.88, 5.83, 5.84, 5.80, 5.89, 5.81, 5.84, 5.83, 5.82, 5.83,
5.81, 5.82, 5.85, 5.81, 5.81, 5.81, 5.84, 5.82, 5.80, 5.84,
5.86, 5.87, 5.82, 5.87, 5.80, 5.81, 5.85, 5.84, 5.83, 5.86,
5.81, 5.81, 5.82, 5.83, 5.85, 5.80, 5.86, 5.82, 5.86, 5.83,
5.80, 5.77, 5.82, 5.85, 5.84, 5.82, 5.85, 5.81, 5.86, 5.79,
5.84, 5.83, 5.80, 5.83, 5.81, 5.83, 5.81, 5.85, 5.83, 5.88,
5.82, 5.87, 5.80, 5.82, 5.83, 5.81, 5.84, 5.79, 5.85, 5.85,
5.84, 5.84, 5.80, 5.82, 5.84, 5.85, 5.86, 5.81, 5.81, 5.85,
5.86, 5.81, 5.81, 5.83, 5.85, 5.85, 5.82, 5.83, 5.86, 5.81 )
LSL=5.65; USL=5.95 
x = data; n=length(x)
#================================================================

# Call Rd3.R from https://github.com/AppliedStat/R-code/blob/master/2020c/Rd3.R
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2020c/Rd3.R")

B=2000;   
#================================================================
Bsample = boot(x, Cpk, R=B, LSL=LSL, USL=USL)  
BCI = boot.ci(Bsample,conf=conf,type= c("norm","basic","perc","bca") )
SCI = SCI.ci(Bsample, conf=conf)

C1 = SCI[2:3]       # SB1
C2 = BCI$norm[2:3]  # SB2
C3 = BCI$perc[4:5]  # PB1
C4 = BCI$basic[4:5] # PB2
C5 = BC.ci(Bsample,conf=conf)[4:5]   # BC1

t.jack=NULL
for ( j in 1L:n ) t.jack = c(t.jack, Cpk(x[-j],LSL=LSL,USL=USL))
BCI1 = boot.ci(Bsample,conf=conf,L=mean(t.jack)-t.jack,type="bca")
C6 = BCI1$bca[4:5]  # BC2

CI = rbind(C1, C2, C3, C4, C5, C6)

colnames(CI) = c("Lower", "Upper")
rownames(CI) = c("SB1", "SB2", "PB1", "PB2", "BC1", "BC2")

cat("\n============================\n")
cat(" Cpk =", Cpk(dat=x,LSL=LSL,USL=USL))
cat("\n Coverage(%)=",conf*100)
cat("\n----------------------------\n\n")

print(CI)


