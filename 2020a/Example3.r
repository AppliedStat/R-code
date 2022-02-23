# -----------------------------------------------------------
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2020a/Rmw3.R")
# -----------------------------------------------------------

# Data from Table 6.4 in Chapter 6 of Montgomery (SQC)
ni = c(5,3,5,5,5, 4,4,5,4,5, 5,5,3,5,3, 5,4,5,5,3, 5,5,5,5,5)
Xbari = c( 74.010, 73.996, 74.008, 74.003, 74.003, 
           73.996, 73.999, 73.997, 74.004, 73.998, 
           73.994, 74.001, 73.994, 73.990, 74.008, 
           73.997, 73.999, 74.007, 73.998, 74.008, 
           74.000, 74.002, 74.002, 74.005, 73.998 ) 
Si    = c( 0.0148, 0.0046, 0.0147, 0.0091, 0.0122, 
           0.0099, 0.0055, 0.0123, 0.0064, 0.0063, 
           0.0029, 0.0042, 0.0100, 0.0153, 0.0087, 
           0.0078, 0.0115, 0.0070, 0.0085, 0.0068, 
           0.0122, 0.0074, 0.0119, 0.0087, 0.0162 ) 
# ------------------------------------------------------------

# ------------------------------------------------------------
nk=3
cat("\n Xbar chart\n")
Xbarchart(Xbari, Si, ni=ni, nk=nk )

cat("\n S chart\n")
Schart( Si, ni, nk=nk )
# ------------------------------------------------------------

# ------------------------------------------------------------
nk=4
cat("\n Xbar chart\n")
Xbarchart(Xbari, Si, ni=ni, nk=nk )

cat("\n S chart\n")
Schart( Si, ni, nk=nk )
# ------------------------------------------------------------

# ------------------------------------------------------------
nk=5
cat("\n Xbar chart\n")
Xbarchart(Xbari, Si, ni=ni, nk=nk )

cat("\n S chart\n")
Schart( Si, ni, nk=nk )
# ------------------------------------------------------------

