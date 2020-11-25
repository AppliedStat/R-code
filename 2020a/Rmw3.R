#--------------------------------------------------------------------
Xbarchart = function(Xbari,Si,ni,nk,CL=c("XB","XA"),FAR=0.002699796){
   CL = match.arg(CL)
   if (CL=="XA") {
      Xbarbar = mean(Xbari)
   }
   else if (CL=="XB") {
      Xbarbar = sum(ni*Xbari) / sum(ni)
   }
   else {
       stop("Choose the Xbarbar: \"XA\" or \"XB\".")
   }
   N = sum(ni)
   m = length(ni)
   c4ni = sqrt(2/(ni-1))*exp(lgamma(ni/2) - lgamma((ni-1)/2))
   one.minus.c4sq = 1-c4ni^2
   S = numeric(4)
   S[1] = sum(Si/c4ni) / m
   S[2] = sum(Si) / sum(c4ni)
   S[3] = sum(c4ni*Si/one.minus.c4sq) / sum(c4ni^2/one.minus.c4sq)
   c4Nm1 = sqrt(2/(N-m))*exp(lgamma((N-m+1)/2) - lgamma((N-m)/2))
   S[4] = sqrt( sum((ni-1)*Si^2) / (N-m) ) / c4Nm1

   z.cut = qnorm(1-FAR/2)
   OUT = array(dim=c(4,3))

   for ( i in 1:4 ) {
        OUT[i,] = Xbarbar + c(-z.cut*S[i]/sqrt(nk),0,z.cut*S[i]/sqrt(nk))
   }
   colnames(OUT) = c("LCL", "CL", "UCL")
   rownames(OUT) = c("A", "B", "C", "D")
   return(OUT)
}
#--------------------------------------------------------------------

#--------------------------------------------------------------------
Schart = function(Si,ni,nk,FAR=0.002699796){ 
   z.cut = qnorm(1-FAR/2)
   c4nk = sqrt(2/(nk-1))*exp(lgamma(nk/2) - lgamma((nk-1)/2))
   c4ni = sqrt(2/(ni-1))*exp(lgamma(ni/2) - lgamma((ni-1)/2))
   N = sum(ni)
   m = length(ni)
   one.minus.c4sq = 1-c4ni^2
   S = numeric(4)
   S[1] = sum(Si/c4ni) / m
   S[2] = sum(Si) / sum(c4ni)
   S[3] = sum(c4ni*Si/one.minus.c4sq) / sum(c4ni^2/one.minus.c4sq)
   c4Nm1 = sqrt(2/(N-m))*exp(lgamma((N-m+1)/2) - lgamma((N-m)/2))
   S[4] = sqrt( sum((ni-1)*Si^2) / (N-m) ) / c4Nm1

   OUT = array(dim=c(4,3))
   for ( i in 1:4 ) {
       CL = c4nk*S[i]
      LCL = max(CL - z.cut*sqrt(1-c4nk^2)*S[i], 0)
      UCL = CL + z.cut*sqrt(1-c4nk^2)*S[i]
      OUT[i,] = c(LCL, CL, UCL)
   }
   colnames(OUT) = c("LCL", "CL", "UCL")
   rownames(OUT) = c("A", "B", "C", "D")
   return(OUT)
}
#--------------------------------------------------------------------

