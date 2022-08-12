#========================================================================
Cpk = function(dat, idx, LSL, USL) {
      if (missing(idx)) idx=c(1:length(dat)) 
      s = sd(dat[idx])
      xbar = mean(dat[idx])
      min( USL-xbar, xbar-LSL ) / (3*s)
}
#------------------------------------------------------------------------
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
