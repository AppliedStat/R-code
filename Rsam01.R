#------------------------------------------------------------------------
within.interval <- function(x, interval){
   stopifnot(length(interval) == 2L)
   interval[1] < x & x < interval[2]
}
#========================================================================

#========================================================================
# Traditional CI 
CI.mean = function(x, conf=0.95) {
   n = length(x)
   t.critical = qt( (1+conf)/2, df=n-1)
   loc = mean(x)
   SE = t.critical * sd(x)/sqrt(n)
   CI = c(loc-SE, loc+SE)
   names(CI) = c("lower", "upper")
   return(CI)
}
#------------------------------------------------------------------------
# Using approx. dist of p-quantile 
CI.qnorm = function(x, p=0.5, conf=0.95) {
   n = length(x)
   z.critical = qnorm( (1+conf)/2 )
   loc = quantile(x, probs=p)
   tmp = sqrt(p*(1-p)/n) * sqrt(2*pi) * sd(x) * exp(0.5*qnorm(p)^2)
   CI = c(loc-z.critical*tmp, loc+z.critical*tmp)
   names(CI) = c("lower", "upper")
   return(CI)
}
#------------------------------------------------------------------------
TI.qnorm = function(x, p=0.5, conf=0.95) {
   n = length(x); xbar = mean(x); s = sd(x)

   delta = -sqrt(n) * qnorm(p) 

   a1 = a2 = (1-conf)/2

   t2 = qt(1-a2, df=n-1, ncp=delta)
   t1 = qt(a1,   df=n-1, ncp=delta)

   CI = c( xbar - t2*s/sqrt(n), xbar - t1*s/sqrt(n) )
   names(CI) = c("lower", "upper")
   return(CI)
}
#------------------------------------------------------------------------
# https://www.itl.nist.gov/div898/handbook/prc/section2/prc263.htm
TI.qnorm.Howe = function(x, p=0.5, conf=0.95) {
   alpha = 1-conf
   n = length(x); xbar = mean(x); s = sd(x)
   qcut = qchisq(alpha, df=n-1)

   u = qnorm( (1+p)/2 ) * sqrt(1+1/n)
   v = sqrt( (n-1)/qcut )
   w = sqrt(1+ (n-3-qcut)/(2*(n+1)^2))

   kcut = u*v*w 

   CI = c(xbar-kcut*s, xbar+kcut*s)
   names(CI) = c("lower", "upper")
   return(CI)
}
#------------------------------------------------------------------------


#=============================================================
# Only one sample. 
# Two-sample case (future work). 
qt.test <-
function (x, alternative = c("two.sided", "less", "greater"),
          xi=0, p=0.5, 
          estimation=c("Default.estimate", "ML.estimate", "MSE.estimate", "UMVU.estimate"), 
          conf.level=0.95 )
{
   alternative <- match.arg(alternative) 
   estimation  <- match.arg(estimation) 
   DNAME = deparse(substitute(x))

   stopifnot(is.numeric(x))
   x = sort(x[complete.cases(x)])
   n =  length(x); df=n-1 

   c4 = sqrt(2/(n-1))*exp(lgamma(n/2) - lgamma((n-1)/2))
   zp = qnorm(p) 
   xbar=mean(x); S=sd(x); sigma.hat = sd(x)*sqrt( (n-1)/n )
   alpha = 1-conf.level

   estimate <- switch(estimation,
          Default.estimate = xbar + zp*S, 
          ML.estimate      = xbar + zp*sigma.hat, 
          MSE.estimate     = xbar + zp*S*c4,
          UMVU.estimate    = xbar + zp*S/c4 )

   names(estimate) = paste( estimation, "of x")

   method = "One sample quantile estimate and noncentral t-test" 

   tstat = (xbar-xi)/(S/sqrt(n))
   noncentrality = -sqrt(n)*zp 

   if (alternative == "less") {
      pval <- pt(tstat, df=df, ncp=noncentrality)
      cint = c(-Inf, xbar-qt(alpha, df=n-1, ncp=noncentrality)*S/sqrt(n))
   }
   else if (alternative == "greater") {
      pval <- pt(tstat, df=df, ncp=noncentrality, lower.tail = FALSE)
      cint = c(xbar-qt(1-alpha, df=n-1, ncp=noncentrality)*S/sqrt(n), Inf)
   }
   else {
      pval1 = pt(tstat, df=df, ncp=noncentrality)
      pval2 = pt(tstat, df=df, ncp=noncentrality, lower.tail = FALSE) 
      pval = 2*min(pval1, pval2)
      t2 = qt(1-alpha/2, df=n-1, ncp=noncentrality)
      t1 = qt(alpha/2,   df=n-1, ncp=noncentrality)
      cint = c( xbar-t2*S/sqrt(n), xbar-t1*S/sqrt(n) )
   }

   names(n)  = "sample size"
   names(xi) = paste(p, "-quantile", sep="")

   names(tstat) = "t" 
   names(noncentrality) = "noncentrality"

   attr(cint, "conf.level") = conf.level
   RVAL = list(statistic=tstat, p.value = pval, parameter= noncentrality, 
               conf.int=cint, sample.size=n, df=n-1,  estimate=estimate, null.value=xi,
               alternative=alternative, method=method, data.name=DNAME)
   class(RVAL) = "htest"
   return(RVAL)
}
#=============================================================

