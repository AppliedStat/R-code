#=====================================================================
# Function for splitting Weibull data and noises using LOO CV idea
weibull.CV = function(data, pvalue=0.05) {
      y = NULL
      test = wp.test(data)
      if ( test$p.value >= pvalue ) return( list(pure=data, outlier=y) )
      while ( test$p.value < pvalue ) {
              test.p.values = numeric(length(data))
              for ( i in 1L:length(data) ) {
                   test.p.values[i] = wp.test(data[-i])$p.value
              }
              j = which.max(test.p.values)
              y = c(y,data[j])
              data = data[-j]
              test = wp.test(data)
      }
      return( list(pure=data, outlier=y) )
}
#---------------------------------------------------------------------
# Function for splitting normal data and noises using boxplot
normal.box = function(data, factor=1.5) {
     iqr = IQR(data)
     quan = quantile(data)
     Q1 = quan[2]
     Q3 = quan[4]
     idx =  (data > Q1-factor*iqr) & (data < Q3+factor*iqr)
     return( list(pure=data[idx], outlier=data[!idx]) )
}
#---------------------------------------------------------------------
# Function for splitting normal data and noises using 3*sigma rule with mean and sd
normal.3sigma.mean = function(data, factor=3) { 
     sigma = sd(data)
     loc   = mean(data)
     idx =  abs(data-loc)/sigma < factor 
     return( list(pure=data[idx], outlier=data[!idx]) )
}
#=====================================================================


#------------------------------
weibull.3sigma.HL = function(data, factor=3, a=0.5) {    # Breakdown point is only 29%
   data = sort(data)
   x = log(data);
   y = log(-log(1-ppoints(x, a=a)))
   DX = outer(x,x,"-"); DY = outer(y,y,"-")
   diag(DX) = NA
   idx = upper.tri(DY)
   alpha.hat = median( (DY[idx]/DX[idx]) )  # slope 
   intercept = median(y-alpha.hat*x) ## This intercept is better than the below. 
   ## intercept = median( apply((outer(x,y,"*")-outer(y,x,"*"))/DX, 2, median, na.rm=TRUE) )
   pred = intercept + alpha.hat*x
   error = y - pred

   sigma = shamos(error)
   idx = abs(error)/sigma < factor
   return( list(pure=data[idx], outlier=data[!idx]) )
}


#------------------------------
weibull.3sigma.med = function(data, factor=3, a=0.5) {    # Breakdown point is only 29%
   data = sort(data)
   x = log(data);
   y = log(-log(1-ppoints(x, a=a)))

   DX = outer(x,x,"-"); DY = outer(y,y,"-")
   diag(DX) = NA
   Sxx = median( apply(DX^2, 2, median, na.rm=TRUE) )   # robust version of Sxx
   Sxy = median( apply(DX*DY,2, median, na.rm=TRUE) )   # robust version of Sxy
   alpha.hat = Sxy / Sxx # slope
   intercept = median(y-alpha.hat*x) ## This intercept is better than the below. 
   ## intercept = median( apply((outer(x,y,"*")-outer(y,x,"*"))/DX, 2, median, na.rm=TRUE) )

   pred = intercept + alpha.hat*x
   error = y - pred

   sigma = mad(error)
   idx = abs(error)/sigma < factor
   return( list(pure=data[idx], outlier=data[!idx]) )
}


#=====================================================================
# Generalized Mean Square Error
gMSE = function(x,y, mux, muy) {
   N = length(x)
   a11 = sum( (x-mux)^2 )
   a22 = sum( (y-muy)^2 )
   a12 = sum( (x-mux)*(y-muy) )
   S = 1/N * matrix( c(a11,a12,a12,a22), nrow=2)
   det(S)
}
#------------------------------
gVAR = function(x,y) {
   det( var(cbind(x,y)) )
}   
#=====================================================================



#=====================================================================
# Shamos estimator
shamos = function (x, constant = 1.048358, na.rm = FALSE, IncludeEqual = FALSE) 
{
    if (na.rm) 
        x <- x[!is.na(x)]
    w1 = outer(x, x, "-")
    w2 = abs(w1[lower.tri(w1, diag = IncludeEqual)])
    constant * median(w2)
}
#=====================================================================
# Hodges-Lehmann estimator
HL = function (x, estimator = c("HL1", "HL2", "HL3"), na.rm = FALSE) 
{
    estimator = match.arg(estimator)
    if (na.rm) 
        x <- x[!is.na(x)]
    xx = outer(x, x, "+")
    HL.estimation = switch(estimator, HL1 = 0.5 * median(xx[lower.tri(xx, 
        diag = FALSE)]), HL2 = 0.5 * median(xx[lower.tri(xx, 
        diag = TRUE)]), HL3 = 0.5 * median(xx))
    return(HL.estimation)
}
#=====================================================================



#=====================================================================
# Random number generation from Gaussian mixture distribution
rweibull.mix = function(n, shapes, scales, prob) {
    npara = length(shapes)
    if (npara != length(scales) ) stop("The numbers of parameters do not match.")
    components = sample (seq_len(npara), prob=prob, size=n, replace=TRUE)
    rweibull(n=n, shape=shapes[components], scale=scales[components])
}
#=====================================================================
# n = 50
# x = rweibull.mix (n, shapes=c(2,5), scales=c(1,5), prob=c(0.9, 0.1)) 







