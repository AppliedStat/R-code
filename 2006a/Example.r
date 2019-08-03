## without log-transform of variance
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2006a/Rsb1.R")
# ===============================================
# DATA INPUT
# -----------------------------------------------
ff <- "CaseStudy.rout"
tau <- 50;

# X data #
u  <- c(-1,0,1)
lu <- length(u)

x1 <- rep(u,lu)
x2 <- rep(u,rep(lu,lu))

X   <- cbind(x1,x2)
l.X <- nrow(X)
nY  <- c(5,5,5, 5,5,5, 5,5,5)

y0.mean <- y0.var <-  numeric(l.X)
yy = array( dim=c(l.X, 5) )
n1 = n2 = n3 = numeric(l.X)
# ===============================================
# Function definition
# -----------------------------------------------
mu     <- function(x) {  50 +  5*(x[1]^2+x[2]^2) }
sigma2 <- function(x) {  10 +  25*(x[1]^2+x[2]^2) }
#
MSE <- function(xt,tau,LM.y, LM.s) {
   ext <- exp(xt)
   x   <- (ext-1)/(ext+1)
   mu  <- predict(LM.y, data.frame(x1=x[1],x2=x[2]) )
   v2  <- predict(LM.s, data.frame(x1=x[1],x2=x[2]) )
   (mu-tau)^2 + v2
}
#
xt <- function(tt)  {ext<-exp(tt); (ext-1)/(ext+1)}
#
xdata <- function(x) data.frame(x1=x[1],x2=x[2])
# ===============================================

mux <- sdx <- numeric(l.X)
for (i in 1:l.X) { mux[i] <- mu(X[i,]); sdx[i] <- sqrt(sigma2(X[i,])) }

# Y data #
set.seed(1)
cat("\n\n*** Data ***\n",file=ff, append=F)   
for (i in 1:l.X) {
     y0 = round(rnorm(nY[i],mux[i],sd=sdx[i]),1)
     yy[i,] = y0 
     tmp = rnorm( 100 ,mux[i],sd=sdx[i])
     n1[i] = sum( tmp < 45 ) 
     n3[i] = sum( tmp > 55 ) 
     n2[i] = 100 - n1[i] - n3[i] 
}
for (i in 1:l.X) {
     y0 = yy[i,]
     X1 = c(y0, rep(-Inf, n1[i]),  rep(45, n2[i]),  rep(55, n3[i]))
     X2 = c(y0, rep(  45, n1[i]),  rep(55, n2[i]),  rep(Inf,n3[i]))
     XX = cbind(X1,X2)
     out = norm.ic.EM (XX, start=c( mean(y0), sqrt(var(y0)))  )

     y0.mean[i] <-  round(out$mu,   3)
     y0.var[i]  <-  round(out$sd^2, 3)
     cat("X =", X[i,], file=ff, append=TRUE)
     cat("\n Y =",  y0, " : mean =", 
         round(mean(y0),2) , "var  =", round(var(y0),2) , file=ff, append=TRUE)
     cat("\n n1 =", n1[i], "n2 =", n2[i], "n3 =", n3[i], file=ff, append=TRUE)
     cat("  :  ", "mean =", y0.mean[i], "var  =", y0.var[i], file=ff, append=TRUE)
     cat("\n\n", file=ff,  append=TRUE)
}
# -----------------------------------------------


# -----------------------------------------------
# Linear Model
# -----------------------------------------------
# -----
# OLS 
# -----
LM.y0mean <-lm(y0.mean ~ (x1+x2)^2 + I(x1^2) + I(x2^2) )
LM.y0var  <-lm(y0.var  ~ (x1+x2)^2 + I(x1^2) + I(x2^2) )

cat("\n\n*** OLS ***\n",file=ff, append=T)   
 
zz <- file(ff, open="at")
sink (file=zz, append=T)
print (summary(LM.y0mean))
print (summary(LM.y0var))

# PREDICT Part
cat("\n\n*** PREDICT ***\n",file=ff, append=T)   
NLM0 <- nlm(MSE,c(0,0),tau=tau,LM.y=LM.y0mean,LM.s=LM.y0var,stepmax=10)
w    <- xt(NLM0$estimate)
mw   <- predict(LM.y0mean,xdata(w))
sw   <- predict(LM.y0var, xdata(w))
cat( " w =", round(w,3), " mean.y =", round(mw,3), " var.y =", round(sw,3), 
     " mean.y - T0 =", round(mw-tau,3), "MSE =", round((mw-tau)^2+sw,3) , "\n",
      file=ff, append=T)


# ===============================================
# CLOSING
# -----------------------------------------------
quit("no")
