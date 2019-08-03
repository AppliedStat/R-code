#=========================================================
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2004/Rk3.R")

## Neet to install cmprsk package
## NB: http://lib.stat.cmu.edu/R/CRAN/web/packages/cmprsk/index.html
##     install.packages("cmprsk")
library(cmprsk)

#===================================================================
# Table IV on Page 16.
#
# Original Data from  W. Nelson (1970)
#      J. of Quality Technology, Vol. 2, No. 3, pg. 126 - 149
# data1 = Group1 (Manual test)
# data2 = Group2 (Manual test)
# data3 = Group3 (Manual test)
#-------------------------------------------------------------------
data1 = c(46,   46, 2145, 2145,  270,   98, 1467,  495, 2145,   16,
          16, 2145,  692, 2145,   52,   98, 2145, 2145,  495, 1961,
        1107, 1937,   12, 1961, 1467,  616, 1453, 1453, 1453,  413,
        1056, 1453, 1453, 1453,  557, 1193)
cause1=c( 6, 3, 0, 0, 6, 6,11,11, 0,10,12, 0,11, 0, 6,11, 0, 0, 6, 0,
        12,11,13, 0,12, 6, 0, 0, 0,11,11, 0, 0, 0,14,11)
#-------------------------------------------------------------------
data2 = c( 311,  571,   73,  571,
        136,  136,  136,  136,  136, 1300, 1300, 1198,  670,  575,
       1300, 1198,  569,  471, 1164, 1164,  608,  608,  608,  608,
        490,  608,  608,  608,  608,   47,  608,  608,  608,  608,
         45,  608,  608,  964,  281,  964,  670, 1164, 1164,  838,
        731,  630,  485,  485,  145,  190,  190)
cause2 = c(11, 0,11,11,
        0, 6, 0, 0, 0, 0, 0, 9,11,11, 0,11, 1,12, 0, 0, 0, 0, 0, 0,
       11,11, 0, 0, 0,11, 0, 0, 0,11, 1, 0, 0,11,12,11,11, 0, 0,11,
        0,11, 0, 0,11, 0, 0)
#-------------------------------------------------------------------
data3 =c( 658,   90,  190,  241,  349,  410,   90,   90,  268,  410,
       410,  485,  508,  631,  631,  631,  635,  658,  731,  739,
       790,  855,  980,  980,  218,  218,  378,  378,  739,  739,
       739,  739,  790,  790,  790,  790,  790,  790,  790,  790,
       790,  790,  790,  980,  980,  980,  980,  980,  600,  600,
       600, 600)
cause3=c( 8, 1, 1, 1, 6,12,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
     11,11,11,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#===================================================================


#-------------------------------------------------------------------
x1 = data1; d1 = cause1
d1[ d1==3 ] = 2; d1[ d1==6 ] = 2; d1[ d1==10] = 2;  d1[ d1==12] = 2;
d1[ d1==13] = 2; d1[ d1==14] = 2; d1[ d1==11] = 1;
#-------------------------------------------------
x2 = data2; d2 = cause2
d2[ d2==1 ] = 2; d2[ d2==6 ] = 2; d2[ d2==9 ] = 2;  d2[ d2==12] = 2;
d2[ d2==11] = 1;
#-------------------------------------------------
x3 = data3; d3 = cause3
d3[ d3==1 ] = 2; d3[ d3==6 ] = 2; d3[ d3==8 ] = 2;  d3[ d3==12] = 2; 
d3[ d3==11] = 1; 
#-------------------------------------------------


#=================================================
pdf( file="ExamplePlot.pdf",  width=7.0, height=8.5)
par( mfrow=c(3,2) )
par(mar=c(5,5,2,2), omi=c(0,0,0,0), cex=0.6, mex=0.8)

#=================================================
# Group 1 
#=================================================
xx  = cuminc (x1,d1)
lam = expo.cr(x1,d1)

cat("\n\n Table I: Parameter Estimates Without Any Restrictions \n")
cat("\n k = 1 \n")
lam

#--------
# Cause 1
#--------
ti = xx[[1]]$time;   P = xx[[1]]$est;  v = xx[[1]]$var 
idx = (P>0)
lower = exp( log(P[idx]) - 1.96*sqrt(v[idx])/P[idx] )
upper = exp( log(P[idx]) + 1.96*sqrt(v[idx])/P[idx] )

plot(range(ti), range(c(lower, upper)) , type="n", xlab="ti", ylab="CIF", sub="(a)", bty="o")
lines(ti, P, lwd=1.3)
lines(ti[idx], lower, lty=2, lwd=1.3)
lines(ti[idx], upper, lty=2, lwd=1.3)

TI = seq( min(ti), max(ti), l=100)
Q  = lam[1]/sum(lam)*( 1- exp(-sum(lam)*TI) )
lines( TI, Q , lwd=1.3)

#--------
# Cause 2
#--------
ti = xx[[2]]$time;   P = xx[[2]]$est;  v = xx[[2]]$var
idx = (P>0)
lower = exp( log(P[idx]) - 1.96*sqrt(v[idx])/P[idx] )
upper = exp( log(P[idx]) + 1.96*sqrt(v[idx])/P[idx] )

plot(range(ti), range(c(lower, upper)) , type="n", xlab="ti", ylab="CIF", sub="(b)", bty="o")
lines(ti, P, lwd=1.3)
lines(ti[idx], lower, lty=2, lwd=1.3)
lines(ti[idx], upper, lty=2, lwd=1.3)

TI = seq( min(ti), max(ti), l=100)
Q  = lam[2]/sum(lam)*( 1- exp(-sum(lam)*TI) )
lines( TI, Q , lwd=1.3)





#=================================================
# Group 2 
#=================================================
xx  = cuminc (x2,d2)
lam = expo.cr(x2,d2)

cat("\n k = 2 \n" )
lam

#--------
# Cause 1
#--------
ti = xx[[1]]$time;   P = xx[[1]]$est;  v = xx[[1]]$var 
idx = (P>0)
lower = exp( log(P[idx]) - 1.96*sqrt(v[idx])/P[idx] )
upper = exp( log(P[idx]) + 1.96*sqrt(v[idx])/P[idx] )

plot(range(ti), range(c(lower, upper)) , type="n", xlab="ti", ylab="CIF", sub="(c)", bty="o")
lines(ti, P, lwd=1.3)
lines(ti[idx], lower, lty=2, lwd=1.3)
lines(ti[idx], upper, lty=2, lwd=1.3)

TI = seq( min(ti), max(ti), l=100)
Q  = lam[1]/sum(lam)*( 1- exp(-sum(lam)*TI) )
lines( TI, Q , lwd=1.3)

#--------
# Cause 2
#--------
ti = xx[[2]]$time;   P = xx[[2]]$est;  v = xx[[2]]$var
idx = (P>0)
lower = exp( log(P[idx]) - 1.96*sqrt(v[idx])/P[idx] )
upper = exp( log(P[idx]) + 1.96*sqrt(v[idx])/P[idx] )

plot(range(ti), range(c(lower, upper)) , type="n", xlab="ti", ylab="CIF", sub="(d)", bty="o")
lines(ti, P, lwd=1.3)
lines(ti[idx], lower, lty=2, lwd=1.3)
lines(ti[idx], upper, lty=2, lwd=1.3)

TI = seq( min(ti), max(ti), l=100)
Q  = lam[2]/sum(lam)*( 1- exp(-sum(lam)*TI) )
lines( TI, Q , lwd=1.3)


#=================================================
# Group 3 
#=================================================
xx  = cuminc (x3,d3)
lam = expo.cr(x3,d3)

cat("\n k = 3 \n" )
lam

#--------
# Cause 1
#--------
ti = xx[[1]]$time;   P = xx[[1]]$est;  v = xx[[1]]$var 
idx = (P>0)
lower = exp( log(P[idx]) - 1.96*sqrt(v[idx])/P[idx] )
upper = exp( log(P[idx]) + 1.96*sqrt(v[idx])/P[idx] )

plot(range(ti), range(c(lower, upper)) , type="n", xlab="ti", ylab="CIF", sub="(e)", bty="o")
lines(ti, P, lwd=1.3)
lines(ti[idx], lower, lty=2, lwd=1.3)
lines(ti[idx], upper, lty=2, lwd=1.3)

TI = seq( min(ti), max(ti), l=100)
Q  = lam[1]/sum(lam)*( 1- exp(-sum(lam)*TI) )
lines( TI, Q , lwd=1.3)

#--------
# Cause 2
#--------
ti = xx[[2]]$time;   P = xx[[2]]$est;  v = xx[[2]]$var
idx = (P>0)
lower = exp( log(P[idx]) - 1.96*sqrt(v[idx])/P[idx] )
upper = exp( log(P[idx]) + 1.96*sqrt(v[idx])/P[idx] )

plot(range(ti), range(c(lower, upper)) , type="n", xlab="ti", ylab="CIF", sub="(f)", bty="o")
lines(ti, P, lwd=1.3)
lines(ti[idx], lower, lty=2, lwd=1.3)
lines(ti[idx], upper, lty=2, lwd=1.3)

TI = seq( min(ti), max(ti), l=100)
Q  = lam[2]/sum(lam)*( 1- exp(-sum(lam)*TI) )
lines( TI, Q , lwd=1.3)



#================================================
# Hypothesis Test 
#================================================
J=2; K=3

#------------------------------------------------
XX = list(x1,x2,x3); dd = list(d1,d2,d3)
 Lam  = expo.cr  (XX,dd)
 Lam00= expo.cr00(XX,dd)   # H01
 Lam10= expo.cr10(XX,dd)   # H02
 Lam01= expo.cr01(XX,dd)   # H03

cat("\n\n Table II: Parameter Estimates Under H_01, H_02, H_03 \n")
cat("\n Under H_01 \n" )
Lam00 

cat("\n Under H_02 \n" )
Lam10 

cat("\n Under H_03 \n" )
Lam01 

 Stat1 = 2 * ( loglike(Lam,XX,dd)-loglike00(Lam00,XX,dd) )
 Stat2 = 2 * ( loglike(Lam,XX,dd)-loglike10(Lam10,XX,dd) )
 Stat3 = 2 * ( loglike(Lam,XX,dd)-loglike01(Lam01,XX,dd) )
#-----------------------------------------------------------
 alpha = 0.05
 CV = c( qchisq(1-alpha,df=J*K-1), qchisq(1-alpha,df=J*K-J), qchisq(1-alpha,df=J*K-K) )

cat("\n\n\n Table III: Test Statistics and Critical Values\n\n")
OUT = cbind( c(Stat1, Stat2, Stat3), c(J*K-1,J*K-J,J*K-K), CV )
colnames(OUT) = c("Test Stat", "df", "Critical Values")
rownames(OUT) = c("S1", "S2", "S3")

OUT 
