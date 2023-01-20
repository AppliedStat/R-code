

source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023c/Wanomaly.R")

# Add weibullness package (https://cran.r-project.org/web/packages/weibullness/)
library(weibullness)

# ===========================================================
# Data set (Reticulum Cell Sarcoma) from
# -----------------------------------------------------------
# 1. Hoel, D. (1972). A representation of mortality data by competing risks. 
#    Biometrics 28, 475–488.
#
# 2. Kalbfleisch, J. and R. Prentice (2002). The statistical analysis of failure time data.
#    Wiley, New York (See Data Set V in Appendix A).
# -----------------------------------------------------------
Data = c(317, 318, 399, 495, 525, 536, 549, 552, 554, 337, 558, 571, 586, 
         594, 596, 605, 612, 621, 628, 631, 636, 643, 647, 648, 649, 661, 
         663, 666, 670, 695, 697, 700, 705, 712, 713, 738, 748, 753)
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




