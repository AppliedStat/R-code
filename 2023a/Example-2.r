
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023a/anomaly.R")
library(MASS)

##
## Short's data set from
## Stigler, S. M. (1977). Do robust estimators work with real data?
## Annals of Statistics, Vol.5, pp.1055-1098.
##
Data <- c( 8.65, 8.35, 8.71, 8.31, 8.36, 8.58, 7.80, 7.71, 8.30, 9.71, 
           8.50, 8.28, 9.87, 8.86, 5.76, 8.44, 8.23) 

# Shapiro test
shapiro.test(Data) # No normal. (p-value is so small)


# ===========================================================================
# Using 3*sigma rule with mean and sd
# ---------------------------------------------------------------------------

# Removing outliers
normal.3sigma.mean(Data)

output = normal.3sigma.mean(Data)

# Data set without outliers
output$pure

# Shapiro test
shapiro.test(output$pure) # Normal data. (p-value is still small)


# ===========================================================================
# Using boxplot
# ---------------------------------------------------------------------------

# Removing outliers
normal.box(Data)

output = normal.box(Data)

# Data set without outliers
output$pure

# Shapiro test
shapiro.test(output$pure) # Normal data. (p-value is large enough)

# ===========================================================================
# Using median and MAD 
# ---------------------------------------------------------------------------

# Removing outliers
normal.median(Data)

# Save the data 
output = normal.median(Data)

# Data set without outliers
output$pure

# Shapiro test
shapiro.test(output$pure) # Normal data. (p-value is large enough)


# ===========================================================================
# Using Huber and RC
# ---------------------------------------------------------------------------

# Removing outliers
normal.huber(Data)

# Save the data 
output = normal.huber(Data)

# Data set without outliers
output$pure

# Shapiro test
shapiro.test(output$pure) # Normal data. (p-value is large enough)


