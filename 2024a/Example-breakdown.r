#
source("https://raw.githubusercontent.com/AppliedStat/R-code/master/2024a/Rmw63.R")

#=======================
# Example (Weibull data)
#=======================
data1 <- c(17.88,  28.92,  33.00,  41.52,  45.12, 45.60, 48.48, 51.84, 51.96,
           54.12,  55.56,  67.80,  68.64,  68.64, 68.88, 84.12, 93.12, 98.64,
          105.12, 105.84, 127.92, 128.04, 173.40)

x = log(sort(data1))
RM.breakdown.point(x,power=1)

RM.breakdown.point(x,power=1/2)

RM.breakdown.point(x,power=1/4)



#=======================
# Example (Birnbaum & Saunders data)
#=======================
data2 = c(3.70,  7.06,  7.16,  7.46,  7.85,  7.97,  8.44,  8.55,  8.58,  8.86,
          8.86,  9.30,  9.60,  9.88,  9.90, 10.00, 10.10, 10.16, 10.18, 10.20, 
         10.55, 10.85, 11.02, 11.02, 11.08, 11.15, 11.20, 11.34, 11.40, 11.99, 
         12,    12,    12.03, 12.22, 12.35, 12.38, 12.52, 12.58, 12.62, 12.69, 
         12.70, 12.90, 12.93, 13.00, 13.10, 13.13, 13.15, 13.30, 13.55, 13.90, 
         14.16, 14.19, 14.20, 14.20, 14.50, 14.52, 14.75, 14.78, 14.81, 14.85, 
         15.02, 15.05, 15.13, 15.22, 15.22, 15.30, 15.40, 15.60, 15.67, 15.78, 
         15.94, 16.02, 16.04, 16.08, 16.30, 16.42, 16.74, 17.30, 17.50, 17.50, 
         17.63, 17.68, 17.81, 17.82, 17.92, 18.20, 18.68, 18.81, 18.90, 18.93, 
         18.95, 19.10, 19.23, 19.40, 19.45, 20.23, 21.00, 21.30, 22.15, 22.68, 24.4)

RM.breakdown.point(data2,power=1/2)

RM.breakdown.point(data2,power=1/4)

RM.breakdown.point(data2,power=1/8)

