#--------------------------------------------------------------------
c4 = function(n) { sqrt(2/(n-1))*exp(lgamma(n/2) - lgamma((n-1)/2)) }
#--------------------------------------------------------------------
Xbarbar = function(Xbari, ni) {
   XA = mean(Xbari)
   XB = sum(ni*Xbari) / sum(ni) 
   return(list(XA=XA, XB=XB))
}
# 
#---------------------------------------
Sbar =  function(Si, ni) {
   N = sum(ni) 
   m = length(ni)
   c4ni = c4(ni)
   one.minus.c4sq = 1-c4ni^2
   SA = sum(Si/c4ni) / m 
   SB = sum(Si) / sum(c4ni) 
   SC = sum(c4ni*Si/one.minus.c4sq) / sum(c4ni^2/one.minus.c4sq)
   SD = sqrt( sum((ni-1)*Si^2) / (N-m) ) / c4(N-m+1)
   Sbar1 = mean(Si)
   Sbar2 = mean(Si) / c4(mean(ni))
   Sbarw = sum(ni*Si) / N  
   OUT=list(SA=SA,SB=SB,SC=SC,SD=SD,Sbar1=Sbar1,Sbar2=Sbar2,Sbarw=Sbarw)
   return(OUT)
}
#--------------------------------------
RE.Sbar = function(ni) {
   N = sum(ni); m = length(ni);
   c4ni = c4(ni)
   varSA = 1 / (m^2) * sum(1/c4ni^2-1)
   varSB = sum(1-c4ni^2)/(sum(c4ni)^2)
   varSC = 1 / sum(c4ni^2/(1-c4ni^2))
   varSD = (1/ c4(N-m+1)^2 - 1)
   varSE = 1/ c4(N)^2 - 1
   OUT = varSE/c(varSA,varSB,varSC,varSD)
   names(OUT) = c("RE(SA)","RE(SB)","RE(SC)","RE(SD)")
   return(OUT)
}
#--------------------------------------
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
   c4ni = c4(ni)
   one.minus.c4sq = 1-c4ni^2
   S = numeric(4)
   S[1] = sum(Si/c4ni) / m
   S[2] = sum(Si) / sum(c4ni)
   S[3] = sum(c4ni*Si/one.minus.c4sq) / sum(c4ni^2/one.minus.c4sq)
   S[4] = sqrt( sum((ni-1)*Si^2) / (N-m) ) / c4(N-m+1)

   z.cut = qnorm(1-FAR/2)
   OUT = array(dim=c(4,3))
   
   for ( i in 1:4 ) {
        OUT[i,] = Xbarbar + c(-z.cut*S[i]/sqrt(nk),0,z.cut*S[i]/sqrt(nk))
   }
   colnames(OUT) = c("LCL", "CL", "UCL")
   rownames(OUT) = c("A", "B", "C", "D")
   return(OUT)
}
#--------------------------------------
Schart = function(Si,ni,nk,FAR=0.002699796){ 
   z.cut = qnorm(1-FAR/2)
   c4nk = c4(nk)
   N = sum(ni)
   m = length(ni)
   c4ni = c4(ni)
   one.minus.c4sq = 1-c4ni^2
   S = numeric(4)
   S[1] = sum(Si/c4ni) / m
   S[2] = sum(Si) / sum(c4ni)
   S[3] = sum(c4ni*Si/one.minus.c4sq) / sum(c4ni^2/one.minus.c4sq)
   S[4] = sqrt( sum((ni-1)*Si^2) / (N-m) ) / c4(N-m+1)

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
#--------------------------------------

