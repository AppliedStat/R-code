
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023a/anomaly.R")
library(MASS)


# Data set from: 
# Wu and Pearn (2008).
# A variables sampling plan based on Cpmk for product acceptance determination.
# European Journal of Operational Research, 184:549-560.
# https://doi.org/10.1016/j.ejor.2006.11.032

Data = c(
0.717,0.698,0.726,0.684,0.727,0.688,0.708,0.703,0.694,0.713,0.730,0.699,0.710,
0.688,0.665,0.704,0.725,0.729,0.716,0.685,0.712,0.716,0.712,0.733,0.709,0.703,
0.730,0.716,0.688,0.688,0.712,0.702,0.726,0.669,0.718,0.714,0.726,0.683,0.713,
0.737,0.740,0.706,0.726,0.688,0.715,0.704,0.724,0.713,0.694,0.742,0.690,0.704,
0.697,0.705,0.707,0.687,0.718,0.718,0.724,0.706,0.687,0.673,0.730,0.732,0.720,
0.688,0.710,0.707,0.706,0.709,0.729,0.792,0.685,0.686,0.722,0.720,0.715,0.727,
0.696)


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
shapiro.test(output$pure) # Normal data. (p-value is large enough)

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
# Using HL amd Shamos
# ---------------------------------------------------------------------------

# Removing outliers
normal.HL(Data)

# Save the data 
output = normal.HL(Data)

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


