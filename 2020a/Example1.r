
##-------------------------------------------------------------------------
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2020a/Rmw3.R")
##-------------------------------------------------------------------------

# Data from  Table~30 in Section 3.31 of \cite{ASTM:2018}
ni = c(50,   50,   100,  25,   25,   50,   100,  50,   50,   50) ## Original
Xbari = c(55.7, 54.6, 52.6, 55.0, 53.4, 55.2, 53.3, 52.3, 53.7, 54.3)
Si    = c(4.35, 4.03, 2.43, 3.56, 3.10, 3.30, 4.18, 4.30, 2.09, 2.67)
##-------------------------------------------------------------------------

cat("\n Xbar chart\n")
Xbarchart(Xbari, Si, ni=ni, nk=25 )

cat("\n S chart\n")
Schart( Si, ni, nk=25 )

cat("\n Xbar chart\n")
Xbarchart(Xbari, Si, ni=ni, nk=50 )

cat("\n S chart\n")
Schart( Si, ni, nk=50 )

cat("\n Xbar chart\n")
Xbarchart(Xbari, Si, ni=ni, nk=100 )

cat("\n S chart\n")
Schart( Si, ni, nk=100 )


