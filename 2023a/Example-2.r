
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2023a/anomaly.R")


# Data set from: 
# Johnson, R. E. and B. H. McFarland (1993).
# Antipsychotic Drug Exposure in a Health Maintenance Organization.
# Medical Care 31(5), 432-444.

Data = c( 43.4, 24,   1.8,  0,   0.1, 170.1,  0.4,  150,  31.5,  5.2,  
          35.7, 27.3, 5,   64.3, 70,     94, 61.9,  9.1,  38.8, 14.8 )

# Shapiro test
shapiro.test(Data) # No normal. (p-value is so small)

# Removing outliers
normal.median(Data)

# Save the data 
output = normal.median(Data)

# Data set without outliers
output$pure

# Shapiro test
shapiro.test(output$pure) # Normal data. (p-value is almost 5%)

#----------------------
# NOTE: The 3-sigma method cannot remove outliers
normal.3sigma.mean(Data)

