# setwd("C:/Users/jm06/Desktop/RCppFunctions")

library(Rcpp)

sourceCpp("CDM1probs.cpp")
sourceCpp("CDM1start1.cpp")
sourceCpp("CDM1start2.cpp")
sourceCpp("CDM1logL.cpp")
sourceCpp("CDM1gr.cpp")
sourceCpp("CDM1hess.cpp")

sourceCpp("CDM5probs.cpp")

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
