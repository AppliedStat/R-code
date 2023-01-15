
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023b/Rcc22.R"

#========================================================================
# Data from Table 3.10 of Lawless (2003). 
# Statistical Models and Methods for Lifetime Data. Wiley, New York.
# Lower ends in [a,b]
X1 = c(
 8, 0, 24, 17, 17, 24, 16, 13, 11, 16,
18, 17, 32, 23, 44, 10, 0, 5, 12, 11,
33, 31, 13, 19, 34, 13, 16, 35, 15, 11,
22, 48, 30, 13, 10, 8, 4, 11, 14, 4,
34, 30, 18, 16, 35, 21, 11 )
# Upper ends in [a,b]
X2 = c(
 12, 22, 31, 27, 23, 30, 24, Inf, 13, 20,
 25, 26, Inf, Inf, 48, 35, 5, 8, 20, Inf,
 40, Inf, 39, 32, Inf, Inf, 24, Inf, 22, 17,
 32, Inf, 34, Inf, 17, 21, 9, Inf, 19, 8,
 Inf, 36, 24, 60, 39, Inf, 20)
XX = cbind(X1,X2)
#========================================================================


# The EM method can find the global maximizer.
weibull.ic.EM (XX, start=c(1,1))



# The Newton-type method fails. 
nlm(neg.loglike.weibull, p=c(1,1), X=XX)

# The Newton-type method is successful with a certain starting value.
nlm(neg.loglike.weibull, p=c(2,30), X=XX)


