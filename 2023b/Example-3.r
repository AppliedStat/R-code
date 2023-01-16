
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023b/Rcc22.R")

# Full observations

#========================================================================
# Data from Example 8.1 on Page 222 of Leemis (2009). Reliability 2nd ed.
Bearings = c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.48, 51.84, 51.96,
             54.12, 55.56, 67.80, 68.64, 68.64, 68.88, 84.12, 93.12, 98.64, 
            105.12,105.84,127.92,128.04,173.40)
XX = cbind(Bearings,Bearings)
#========================================================================


# The EM method can find the global maximizer.
weibull.ic.EM (XX, start=c(1,1))

weibull.ic.EM (XX, start=c(5,5))


# The Newton-type method fails . 
# Also, the estimates depend on a starting value.
nlm(neg.loglike.weibull, p=c(1,1), X=XX)

nlm(neg.loglike.weibull, p=c(5,5), X=XX)


