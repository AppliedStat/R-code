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
gVAR = function(x,y) { det( var(cbind(x,y)) ) }   
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


