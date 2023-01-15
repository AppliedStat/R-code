
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023b/Rcc22.R")

#========================================================================
# Data from Table 3.10 (Radiotherapy) of Lawless (2003). 
# Statistical Models and Methods for Lifetime Data. Wiley, New York.
# Lower ends in [a,b]
X1 = c(
45, 25, 37, 6, 46, 0, 0, 26, 18, 46,
46, 24, 46, 27, 36, 7, 36, 5, 17, 46,
19, 7, 36, 17, 37, 37, 24, 0, 40, 32,
 4, 17, 33, 15, 46, 19, 11, 11, 37, 22,
38, 34, 46, 5, 36, 46, 14) 
X2 = c(
Inf,  37, Inf,  10, Inf,   5,   7, 40, Inf, Inf,
Inf, Inf, Inf,  34, Inf,  16,  44, 11, Inf, Inf,
 35,  14,  48,  25,  44, Inf, Inf,  8, Inf, Inf,
 11,  25, Inf, Inf, Inf,  26,  15, 18, Inf, Inf,
Inf, Inf, Inf,  12, Inf, Inf,  17)
XX = cbind(X1,X2)
#========================================================================


# The EM method can find the global maximizer.
weibull.ic.EM (XX, start=c(1,1))



# The Newton-type method fails. 
nlm(neg.loglike.weibull, p=c(1,1), X=XX)

# The Newton-type method is successful with a certain starting value.
nlm(neg.loglike.weibull, p=c(2,50), X=XX)


