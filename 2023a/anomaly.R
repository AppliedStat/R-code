#=====================================================================
# Function for splitting normal data and noises using LOO CV idea
normal.CV = function(data, pvalue=0.05) { 
      y = NULL
      test = shapiro.test(data)
      if ( test$p.value >= pvalue ) return( list(pure=data,outlier=y) )
      test.p.values = numeric(length(data))
      while ( test$p.value < pvalue ) {
              test.p.values = numeric(length(data))
	      for ( i in 1L:length(data) ) {
		   test.p.values[i] = shapiro.test(data[-i])$p.value
	      }
	      j = which.max(test.p.values)
	      y = c(y,data[j])
	      data = data[-j]
	      test = shapiro.test(data)
      }
      return( list(pure=data, outlier=y) )
}
#---------------------------------------------------------------------
# Function for splitting normal data and noises using 3*sigma rule with mean and sd
normal.3sigma.mean = function(data, factor=3) { 
     sigma = sd(data)
     loc   = mean(data)
     idx =  abs(data-loc)/sigma < factor 
     return( list(pure=data[idx], outlier=data[!idx]) )
}
#---------------------------------------------------------------------
# Function for splitting normal data and noises using 3*sigma rule with HL and Shamos estimators
normal.HL = function(data, factor=3) { 
     sigma = shamos(data)
     loc   = HL(data)
     idx =  abs(data-loc)/sigma < factor 
     return( list(pure=data[idx], outlier=data[!idx]) )
}
#---------------------------------------------------------------------
# Function for splitting normal data and noises using 3*sigma rule with median and MAD estimators
normal.median = function(data, factor=3) {     
     sigma = mad(data)
     loc   = median(data)
     idx =  abs(data-loc)/sigma < factor
     return( list(pure=data[idx], outlier=data[!idx]) )
}
#---------------------------------------------------------------------
# Function for splitting normal data and noises using 3*sigma rule with Huber and Qn estimators
normal.huber = function(data, factor=3, k=1.5) { 
     sigma = RC(data)
     loc   = huber(data,k=k)$mu
     idx =  abs(data-loc)/sigma < factor
     return( list(pure=data[idx], outlier=data[!idx]) )
}
#=====================================================================


#=====================================================================
# Function for splitting normal data and noises using boxplot
normal.box = function(data, factor=1.5) { 
     iqr = IQR(data)
     quan = quantile(data)
     Q1 = quan[2]
     Q3 = quan[4]
     idx =  (data > Q1-factor*iqr) & (data < Q3+factor*iqr)
     return( list(pure=data[idx], outlier=data[!idx]) )
}
#=====================================================================


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
# Qn or RC estimator
RC <- function(x,d=2.2191444659850758647, IncludeEqual=FALSE) {
   w1 = outer(x,x, "-")
   w2 = abs( w1[lower.tri(w1,diag=IncludeEqual)]  )
   d*quantile(w2,1/4,type=5)
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
rnorm.mix = function(n, means, sds, prob) {
    npara = length(means)
    if (npara != length(sds) ) stop("The numbers of parameters do not match.")
    components = sample (1L:npara, prob=prob, size=n, replace=TRUE)
    rnorm(n=n, mean=means[components], sd=sds[components])
}
#=====================================================================







