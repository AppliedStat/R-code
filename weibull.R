#
# MLE of Weibull (../PP/pp5-EM/pgm)
#
# See Appendix D (Page 338 of Leemis). '
# Refer to Farnum & Booth (1997) Uniqueness of MLE of Weibull. IEEE Rel. 46, 523-525.
# 
weibull.MLE <- function(x, interval) {
   if ( any (x <= 0) ) stop("The data should be positive")
   if (missing(interval)) {
      meanlog = mean(log(x))
      lower = 1 / ( log(max(x)) - meanlog )
      upper = sum( (x^lower)*log(x) ) / sum( x^lower ) - meanlog
      interval = c(lower,1/upper)
   }
   EE = function(alpha,x) {
      xalpha = x^alpha
      sum(log(x)*(xalpha)) / sum(xalpha) - 1/alpha - mean(log(x))
   }
   tmp = uniroot(EE, interval=interval, x=x)
   alpha = tmp$root
   list ( alpha=alpha, lam=1/mean(x^alpha) )
}
# 
# Note: lam = 1/beta (the pdf of the haoudout is used)
#       if beta = scale (R:xweibull), lam=1/(beta^alpha)
#       if lambda is defined as in Lemmis book, lam = lambda^alpha and alpha=kappa



# For Elsayed definition of Weibull pdf
#  Convert (alpha, lam) to (alpha, beta ) parameters in Leemis book
Reparametrization = function(alpha, lam) {
    lambda = lam^(1/alpha)
    beta   = 1/lambda 
    list(alpha=alpha, beta=beta)
}


