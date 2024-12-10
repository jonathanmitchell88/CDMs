# setwd("C:/Users/jm06/Desktop/RCppFunctions")

library(Rcpp)

# Rcpp functions for computing probabilities, (-2log) likelihoods, gradients, Hessians and starting guesses for maximum likelihood estimates of parameters.
# Hessian for CDM 5 to come and functions for CDMs 2 - 4 to come.
# Functions have been tested by comparing to values computed using Mathematica.

sourceCpp("CDM1probs.cpp")
sourceCpp("CDM1start1.cpp")
sourceCpp("CDM1start2.cpp")
sourceCpp("CDM1logL.cpp")
sourceCpp("CDM1gr.cpp")
sourceCpp("CDM1hess.cpp")

sourceCpp("CDM5probs.cpp")
sourceCpp("CDM5start.cpp")
sourceCpp("CDM5logL.cpp")
sourceCpp("CDM5gr.cpp")

#

n <- 1e6

l <- 1e4
m <- matrix(0,l,21)
iters1 <- rep(0,l)
iters2 <- rep(0,l)
iters3 <- rep(0,l)

start.time <- Sys.time()

for (i in 1:l) {
  pars <- c(runif(1,-1,1),runif(5))
  probs <- CDM1probsRcpp(pars,'()')
  data <- rmultinom(1,n,probs)
  
  start11 <- CDM1start1Rcpp(data,'()',1e-2)
  start21 <- CDM1start2Rcpp(data,'()',1e-2)
  like11 <- CDM1logLRcpp(start11,data,'()')
  like21 <- CDM1logLRcpp(start21,data,'()')
  
  start12 <- CDM1start1Rcpp(data,'(23)',1e-2)
  start22 <- CDM1start2Rcpp(data,'(23)',1e-2)
  like12 <- CDM1logLRcpp(start12,data,'(23)')
  like22 <- CDM1logLRcpp(start22,data,'(23)')
  
  start13 <- CDM1start1Rcpp(data,'(234)',1e-2)
  start23 <- CDM1start2Rcpp(data,'(234)',1e-2)
  like13 <- CDM1logLRcpp(start13,data,'(234)')
  like23 <- CDM1logLRcpp(start23,data,'(234)')
  
  if (like11<=like21) {
    summary1 <- nlminb(start=start11,objective=CDM1logLRcpp,gradient=CDM1grRcpp,hessian=CDM1hessRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='()')
  } else {
    summary1 <- nlminb(start=start21,objective=CDM1logLRcpp,gradient=CDM1grRcpp,hessian=CDM1hessRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='()')
  }
  
  if (like12<=like22) {
    summary2 <- nlminb(start=start12,objective=CDM1logLRcpp,gradient=CDM1grRcpp,hessian=CDM1hessRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(23)')
  } else {
    summary2 <- nlminb(start=start22,objective=CDM1logLRcpp,gradient=CDM1grRcpp,hessian=CDM1hessRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(23)')
  }
  
  if (like13<=like23) {
    summary3 <- nlminb(start=start13,objective=CDM1logLRcpp,gradient=CDM1grRcpp,hessian=CDM1hessRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(234)')
  } else {
    summary3 <- nlminb(start=start23,objective=CDM1logLRcpp,gradient=CDM1grRcpp,hessian=CDM1hessRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(234)')
  }
  
  m[i,1:6] <- summary1$par
  m[i,7] <- summary1$objective*n
  
  m[i,8:13] <- summary2$par
  m[i,14] <- summary2$objective*n
  
  m[i,15:20] <- summary3$par
  m[i,21] <- summary3$objective*n
  
  iters1[i] <- summary1$iterations
  iters2[i] <- summary2$iterations
  iters3[i] <- summary3$iterations
}

end.time <- Sys.time()
time.taken <- end.time-start.time

time.taken

c(min(iters1),mean(iters1),max(iters1))
c(min(iters2),mean(iters2),max(iters2))
c(min(iters3),mean(iters3),max(iters3))

count <- 0

for (i in 1:l) {
  if (which.min(m[i,c(7,14,21)])==1) {
    count <- count+1
  }
}

count/l

#

n <- 1e6

l <- 1e4
m <- matrix(0,l,66)
iters1 <- rep(0,l)
iters2 <- rep(0,l)
iters3 <- rep(0,l)
iters4 <- rep(0,l)
iters5 <- rep(0,l)
iters6 <- rep(0,l)

start.time <- Sys.time()

for (i in 1:l) {
  pars <- c(runif(1,-1,1),runif(9))
  probs <- CDM5probsRcpp(pars,'()')
  data <- rmultinom(1,n,probs)
  
  summary1 <- nlminb(start=CDM5startRcpp(data,'()',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='()')
  summary2 <- nlminb(start=CDM5startRcpp(data,'(34)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(34)')
  summary3 <- nlminb(start=CDM5startRcpp(data,'(23)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(23)')
  summary4 <- nlminb(start=CDM5startRcpp(data,'(243)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(243)')
  summary5 <- nlminb(start=CDM5startRcpp(data,'(234)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(234)')
  summary6 <- nlminb(start=CDM5startRcpp(data,'(24)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(24)')
  
  m[i,1:10] <- summary1$par
  m[i,11] <- summary1$objective*n
  
  m[i,12:21] <- summary2$par
  m[i,22] <- summary2$objective*n
  
  m[i,23:32] <- summary3$par
  m[i,33] <- summary3$objective*n
  
  m[i,34:43] <- summary4$par
  m[i,44] <- summary4$objective*n
  
  m[i,45:54] <- summary5$par
  m[i,55] <- summary5$objective*n
  
  m[i,56:65] <- summary6$par
  m[i,66] <- summary6$objective*n
  
  iters1[i] <- summary1$iterations
  iters2[i] <- summary2$iterations
  iters3[i] <- summary3$iterations
  iters4[i] <- summary4$iterations
  iters5[i] <- summary5$iterations
  iters6[i] <- summary6$iterations
}

end.time <- Sys.time()
time.taken <- end.time-start.time

time.taken

c(min(iters1),mean(iters1),max(iters1))
c(min(iters2),mean(iters2),max(iters2))
c(min(iters3),mean(iters3),max(iters3))
c(min(iters4),mean(iters4),max(iters4))
c(min(iters5),mean(iters5),max(iters5))
c(min(iters6),mean(iters6),max(iters6))

count <- 0

for (i in 1:l) {
  if (which.min(m[i,c(11,22,33,44,55,66)])==1) {
    count <- count+1
  }
}

count/l
