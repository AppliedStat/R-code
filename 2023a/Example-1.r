
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023a/anomaly.R")
library(MASS)

##
##- Newcomb's light speed data set: 
## Stigler, S. M. (1977). Do robust estimators work with real data?
## Annals of Statistics, Vol.5, pp.1055-1098.
##
Data <- c(28,-44,29,30,26,27,22,23,33,16,24,29,24,40,21,31,34,-2,25,19,24,28,
          37,32,20,25,25,36,36,21,28,26,32,28,26,30,36,29,30,22,36,27,26,28,
          29,23,31,32,24,27,27,27,32,25,28,27,26,24,32,29,28,33,39,25,16,23)

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


