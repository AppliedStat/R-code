#========================================================================
Cpk.true = function(mu, sigma, LSL, USL) {pmin(USL-mu,mu-LSL)/(3*sigma)}
#------------------------------------------------------------------------
Cpk = function(dat, idx, LSL, USL) {
      if (missing(idx)) idx=c(1:length(dat)) 
      s = sd(dat[idx])
      xbar = mean(dat[idx])
      min( USL-xbar, xbar-LSL ) / (3*s)
}
#------------------------------------------------------------------------
rCpk = function(dat, idx, LSL, USL) {
      if (missing(idx)) idx=c(1:length(dat))
      s = mad(dat[idx])          ## MAD is used.
      xbar = median(dat[idx])    ## median is used.
      min( USL-xbar, xbar-LSL ) / (3*s)
}
#------------------------------------------------------------------------
# Expectation of Cpk.hat. See Lin/Pearn (2003). Dist. of Estimate of Cpk.
ECpk = function(n, mu, sigma, LSL, USL) {
   m = (USL+LSL)/2; d = (USL-LSL)/2 
  Cp = (USL-LSL)/6/sigma ; Ca = 1 - abs(mu-m)/d 
  bnminus1 = sqrt(2/(n-1))*exp(lgamma((n-1)/2)-lgamma((n-2)/2))
  lam = 9*n*(Cp-Cp*Ca)^2 
  tmp = abs(mu-m)/sigma*(1-2*pnorm(-sqrt(lam)))
  (d/sigma - sqrt(2/n/pi)*exp(-lam/2) - tmp)/(3*bnminus1)
}
#------------------------------------------------------------------------
# CI1 (approx CI for Cpk based on N(0,1): Mathews (2010)
#      See also Montgomery EQ (8.21) on  Page 371.
CI1 = function(dat, LSL, USL, conf) {
    n  = length(dat)
    Cpk.hat = Cpk(dat=dat, LSL=LSL, USL=USL)
    delta = sqrt( 1/n * (1/9/Cpk.hat^2+0.5) )
    z.cut = qnorm( (1+conf)/2 )
    c( Cpk.hat*(1-z.cut*delta),  Cpk.hat*(1+z.cut*delta) )
}
# CI2 (robustified version)
rCI1 = function(dat, LSL, USL, conf) {
    n  = length(dat)
    rCpk.hat = rCpk(dat=dat, LSL=LSL, USL=USL)
    delta = sqrt( 1/n * (1/9/rCpk.hat^2+0.5) )
    z.cut = qnorm( (1+conf)/2  )
    c( rCpk.hat*(1-z.cut*delta), rCpk.hat*(1+z.cut*delta) )
}

#==========================================================================
Cpc.true = function(sigma, LSL, USL) {(USL-LSL)/(6*sigma)} # Under Normal 
Cpc = function(dat, idx, LSL, USL) {
      if (missing(idx)) idx=c(1:length(dat))
      m = (LSL+USL)/2
      (USL-LSL) / ( 7.5198848238930011689*mean(abs(dat[idx]-m)) )
}
#==========================================================================



#==========================================================================
#==========================================================================
#==========================================================================
SCI.ci = function (boot.out = NULL, conf = 0.95, index = 1, var.t0 = NULL, 
    t0 = NULL, t = NULL, L = NULL, h = function(t) t, hdot = function(t) 1, 
    hinv = function(t) t) 
{
    if (is.null(t0)) {
        if (!is.null(boot.out)) 
            t0 <- boot.out$t0[index]
        else stop("bootstrap output object or 't0' required")
    }
    if (!is.null(boot.out) && is.null(t)) 
        t <- boot.out$t[, index]
    if (!is.null(t)) {
        fins <- seq_along(t)[is.finite(t)]
        t <- h(t[fins])
    }
    if (is.null(var.t0)) var.t0 <- var(t) 
    else var.t0 <- var.t0 * hdot(t0)^2
    t0 <- h(t0)
    bias = 0
    merr <- sqrt(var.t0) * qnorm((1 + conf)/2)
    out <- cbind(conf, hinv(t0-bias-merr), hinv(t0 -bias+merr))
    out
}
#------------------------------------------------------------------------
# See boot:::bca.ci
BC.ci = function (boot.out, conf=0.95, index=1, t0=NULL, t=NULL, L=NULL,
    h = function(t) t, hdot = function(t) 1, hinv = function(t) t, ...) 
{
    t.o <- t
    if (is.null(t) || is.null(t0)) {
        t <- boot.out$t[, index]
        t0 <- boot.out$t0[index]
    }
    t <- t[is.finite(t)]
    w <- qnorm(sum(t < t0)/length(t))
    if (!is.finite(w)) 
        stop("estimated adjustment 'w' is infinite")
    alpha <- (1 + c(-conf, conf))/2
    zalpha <- qnorm(alpha)
    a = 0  # 
    adj.alpha <- pnorm(w + (w + zalpha)/(1 - a * (w + zalpha)))
    qq <- boot:::norm.inter(t, adj.alpha)
    cbind(conf, matrix(qq[,1L], ncol=2L), matrix(hinv(h(qq[,2L])), ncol=2L))
}
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
within.interval <- function(x, interval){
   stopifnot(length(interval) == 2L)
   interval[1] < x & x < interval[2]
}
#===========================================================================
Boot.CI1 = function (x, statistic, conf=conf, R=1000, LSL, USL) 
{
     Bsample = boot(x, statistic=statistic, R=R, LSL=LSL, USL=USL)
   # BCI: norm, basic, perc, bca 
     BCI = boot.ci(Bsample,conf=conf,type= c("norm","basic","perc","bca") )

   # SCI: Franklin/Wasserman
     SCI = SCI.ci(Bsample, conf=conf)

     C1 = SCI[2:3]       # SB1
     C2 = BCI$norm[2:3]  # SB2
     C3 = BCI$perc[4:5]  # PB1
     C4 = BCI$basic[4:5] # PB2
     C5 = BC.ci(Bsample, conf=conf)[4:5]   # BC1
  
     t.jack=NULL
     for ( j in 1L:length(x) )  {
         t.jack = c(t.jack, statistic(x[-j],LSL=LSL,USL=USL))
     }
     BCI1 = boot.ci(Bsample,conf=conf,L=mean(t.jack)-t.jack,type="bca")
     C6 = BCI1$bca[4:5]  # BC2
     CI = rbind(C1, C2, C3, C4, C5, C6)
     CI = cbind(CI, CI[,2]-CI[,1] )
     colnames(CI) = c("lower", "upper", "width")
     rownames(CI) = c("SB1", "SB2", "PB1", "PB2", "BC1", "BC2" )
     return(CI)
}
#---------------------------------------------------------------------------
Boot.CI2 = function (x, statistic, conf=conf, R=1000, LSL, USL) 
{
     Bsample = boot(x, statistic=statistic, R=R, LSL=LSL, USL=USL)
   # BCI: norm, basic, perc, bca 
     BCI = boot.ci(Bsample,conf=conf,type= c("norm","basic","perc","bca") )

   # SCI: Franklin/Wasserman
     SCI = SCI.ci(Bsample, conf=conf)

     C1 = SCI[2:3]       # SB1
     C2 = BCI$norm[2:3]  # SB2
     C3 = BCI$perc[4:5]  # PB1
     C4 = BCI$basic[4:5] # PB2
     C5 = BC.ci(Bsample, conf=conf)[4:5]   # BC1

     CI = rbind(C1, C2, C3, C4, C5)
     CI = cbind(CI, CI[,2]-CI[,1] )
     colnames(CI) = c("lower", "upper", "width")
     rownames(CI) = c("SB1", "SB2", "PB1", "PB2", "BC1" )
     return(CI)
}

#---------------------------------------------------------------------------
