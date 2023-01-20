

source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023c/Wanomaly.R")

# Add weibullness package (https://cran.r-project.org/web/packages/weibullness/)
library(weibullness)

# ===========================================================
# Data set (Other causes) from
# -----------------------------------------------------------
# 1. Hoel, D. (1972). A representation of mortality data by competing risks. 
#    Biometrics 28, 475–488.
#
# 2. Kalbfleisch, J. and R. Prentice (2002). The statistical analysis of failure time data.
#    Wiley, New York (See Data Set V in Appendix A).
# -----------------------------------------------------------
Data = c( 40,  42,  51,  62, 163, 179, 206, 222, 228, 252, 249, 282, 324,
         333, 341, 366, 385, 407, 420, 431, 441, 461, 462, 482, 517, 517,
         524, 564, 567, 586, 619, 620, 621, 622, 647, 651, 686, 761, 763)
# -----------------------------------------------------------

# Weibullness test with the original data
wp.test(Data)

# ---------------------------------------
# Remove outliers using CV method
weibull.CV(Data)
output = weibull.CV(Data)

# Weibullness test again with pure data set
wp.test(output$pure)


# ---------------------------------------
# Remove outliers using 3 sigma method
normal.3sigma.mean(Data)
output = normal.3sigma.mean(Data)

# Weibullness test again with pure data set
wp.test(output$pure)


# ---------------------------------------
# Remove outliers using box-plot method
normal.box(Data)
output = normal.box(Data)

# Weibullness test again with pure data set
wp.test(output$pure)




