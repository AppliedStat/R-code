
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2006b/Rpp5.R")

 X = c(1.9, 2.1, 3.2, 1.1, 2.1, 1.0, 2.0, 6.1, 3)     # lifetime observation
 M = list(1, 1, 1, 2, 2, 3, 3, 0, c(1,2,3))           # failure modes
 d =    c(1, 1, 1, 2, 2, 3, 3, 0, -1)                 # same as M


#  ---------------------------
#  We assume that there is full masking in M. 
   expo.cm.EM(X,M)      # Exponential Model
   norm.cm.EM(X,M)      # Normal Model
   norm.cm.EM(log(X),M) # Lognormal Model
   weibull.cm.EM(X,M)   # Weibull Model
   wald.cm.EM(X,M)      # Wald (inverse Gaussian) Model (need to be improved)

#  ---------------------------
#  We assume that there is partial masking in M. 
   M = list(1, 1, 0, c(2,3), 2, 3, 3, c(1,2), c(1,2,3)) 

   expo.cm.EM(X,M)
   norm.cm.EM(X,M)      # Normal Model
   norm.cm.EM(log(X),M) # Lognormal Model
   weibull.cm.EM(X,M)   # Weibull Model
   wald.cm.EM(X,M)      # Wald (inverse Gaussian) Model (need to be improved)

