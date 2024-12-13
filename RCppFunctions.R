# setwd("C:/Users/jm06/Desktop/RCppFunctions")

library(Rcpp)

# Rcpp functions for computing probabilities, (-2log) likelihoods, gradients, Hessians and starting guesses for maximum likelihood estimates of parameters.
# Functions have been tested by comparing to values computed using Mathematica.

# Estimated times for computing all likelihoods on Jonathan's computer:

# 1 core

# 9-taxon CDM: 0.369s
# 10-taxon CDM: 0.554s
# 20-taxon CDM: 6.39s
# 50-taxon CDM: 2.03 mins
# 100-taxon CDM: 17.2 mins
# 200-taxon CDM: 2.37 hours
# 500-taxon CDM: 1.57 days
# 1000-taxon CDM: 12.6 days

# 8 cores

# 9-taxon CDM: 0.0462s
# 10-taxon CDM: 0.0692s
# 20-taxon CDM: 0.799s
# 50-taxon CDM: 15.2s
# 100-taxon CDM: 2.15 mins
# 200-taxon CDM: 17.8 mins
# 500-taxon CDM: 4.71 hours
# 1000-taxon CDM: 1.58 days

# 100 cores

# 9-taxon CDM: 0.000369s
# 10-taxon CDM: 0.00554s
# 20-taxon CDM: 0.0639s
# 50-taxon CDM: 1.21s
# 100-taxon CDM: 10.3s
# 200-taxon CDM: 1.42 mins
# 500-taxon CDM: 22.6 mins
# 1000-taxon CDM: 3.03 hours

# 9 taxa, 8 cores, 20000 genes

# 15.4 mins

sourceCpp("CDM1probs.cpp")
sourceCpp("CDM1start1.cpp")
sourceCpp("CDM1start2.cpp")
sourceCpp("CDM1logL.cpp")
sourceCpp("CDM1gr.cpp")
sourceCpp("CDM1hess.cpp")

sourceCpp("CDM2probs.cpp")
sourceCpp("CDM2start.cpp")
sourceCpp("CDM2logL.cpp")
sourceCpp("CDM2gr.cpp")

sourceCpp("CDM3probs.cpp")
sourceCpp("CDM3start.cpp")
sourceCpp("CDM3logL.cpp")
sourceCpp("CDM3gr.cpp")

sourceCpp("CDM4probs.cpp")
sourceCpp("CDM4start.cpp")
sourceCpp("CDM4logL.cpp")
sourceCpp("CDM4gr.cpp")

sourceCpp("CDM5probs.cpp")
sourceCpp("CDM5start.cpp")
sourceCpp("CDM5logL.cpp")
sourceCpp("CDM5gr.cpp")

# CDM 1

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
    summary1 <- nlminb(start=start11,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='()')
  } else {
    summary1 <- nlminb(start=start21,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='()')
  }
  
  if (like12<=like22) {
    summary2 <- nlminb(start=start12,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(23)')
  } else {
    summary2 <- nlminb(start=start22,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(23)')
  }
  
  if (like13<=like23) {
    summary3 <- nlminb(start=start13,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(234)')
  } else {
    summary3 <- nlminb(start=start23,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(234)')
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

# CDM 2

n <- 1e6

l <- 1e4
m <- matrix(0,l,48)
iters1 <- rep(0,l)
iters2 <- rep(0,l)
iters3 <- rep(0,l)
iters4 <- rep(0,l)
iters5 <- rep(0,l)
iters6 <- rep(0,l)

start.time <- Sys.time()

for (i in 1:l) {
  pars <- c(runif(1,-1,1),runif(6))
  probs <- CDM2probsRcpp(pars,'()')
  data <- rmultinom(1,n,probs)
  
  summary1 <- nlminb(start=CDM2startRcpp(data,'()',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='()')
  summary2 <- nlminb(start=CDM2startRcpp(data,'(34)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(34)')
  summary3 <- nlminb(start=CDM2startRcpp(data,'(23)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(23)')
  summary4 <- nlminb(start=CDM2startRcpp(data,'(243)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(243)')
  summary5 <- nlminb(start=CDM2startRcpp(data,'(234)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(234)')
  summary6 <- nlminb(start=CDM2startRcpp(data,'(24)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(24)')
  
  m[i,1:7] <- summary1$par
  m[i,8] <- summary1$objective*n
  
  m[i,9:15] <- summary2$par
  m[i,16] <- summary2$objective*n
  
  m[i,17:23] <- summary3$par
  m[i,24] <- summary3$objective*n
  
  m[i,25:31] <- summary4$par
  m[i,32] <- summary4$objective*n
  
  m[i,33:39] <- summary5$par
  m[i,40] <- summary5$objective*n
  
  m[i,41:47] <- summary6$par
  m[i,48] <- summary6$objective*n
  
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
  if (which.min(m[i,c(8,16,24,32,40,48)])==1) {
    count <- count+1
  }
}

count/l

# CDM 3

n <- 1e6

l <- 1e4
m <- matrix(0,l,54)
iters1 <- rep(0,l)
iters2 <- rep(0,l)
iters3 <- rep(0,l)
iters4 <- rep(0,l)
iters5 <- rep(0,l)
iters6 <- rep(0,l)

start.time <- Sys.time()

for (i in 1:l) {
  pars <- c(runif(1,-1,1),runif(9))
  probs <- CDM3probsRcpp(pars,'()')
  data <- rmultinom(1,n,probs)
  
  summary1 <- nlminb(start=CDM3startRcpp(data,'()',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='()')
  summary2 <- nlminb(start=CDM3startRcpp(data,'(34)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(34)')
  summary3 <- nlminb(start=CDM3startRcpp(data,'(23)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(23)')
  summary4 <- nlminb(start=CDM3startRcpp(data,'(243)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(243)')
  summary5 <- nlminb(start=CDM3startRcpp(data,'(234)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(234)')
  summary6 <- nlminb(start=CDM3startRcpp(data,'(24)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(24)')
  
  m[i,1:8] <- summary1$par
  m[i,9] <- summary1$objective*n
  
  m[i,10:17] <- summary2$par
  m[i,18] <- summary2$objective*n
  
  m[i,19:26] <- summary3$par
  m[i,27] <- summary3$objective*n
  
  m[i,28:35] <- summary4$par
  m[i,36] <- summary4$objective*n
  
  m[i,37:44] <- summary5$par
  m[i,45] <- summary5$objective*n
  
  m[i,46:53] <- summary6$par
  m[i,54] <- summary6$objective*n
  
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
  if (which.min(m[i,c(9,18,27,36,45,54)])==1) {
    count <- count+1
  }
}

count/l

# CDM 4

n <- 1e6

l <- 1e4
m <- matrix(0,l,60)
iters1 <- rep(0,l)
iters2 <- rep(0,l)
iters3 <- rep(0,l)
iters4 <- rep(0,l)
iters5 <- rep(0,l)
iters6 <- rep(0,l)

start.time <- Sys.time()

for (i in 1:l) {
  pars <- c(runif(1,-1,1),runif(9))
  probs <- CDM4probsRcpp(pars,'()')
  data <- rmultinom(1,n,probs)
  
  summary1 <- nlminb(start=CDM4startRcpp(data,'()',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='()')
  summary2 <- nlminb(start=CDM4startRcpp(data,'(34)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(34)')
  summary3 <- nlminb(start=CDM4startRcpp(data,'(23)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(23)')
  summary4 <- nlminb(start=CDM4startRcpp(data,'(243)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(243)')
  summary5 <- nlminb(start=CDM4startRcpp(data,'(234)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(234)')
  summary6 <- nlminb(start=CDM4startRcpp(data,'(24)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(24)')
  
  m[i,1:9] <- summary1$par
  m[i,10] <- summary1$objective*n
  
  m[i,11:19] <- summary2$par
  m[i,20] <- summary2$objective*n
  
  m[i,21:29] <- summary3$par
  m[i,30] <- summary3$objective*n
  
  m[i,31:39] <- summary4$par
  m[i,40] <- summary4$objective*n
  
  m[i,41:49] <- summary5$par
  m[i,50] <- summary5$objective*n
  
  m[i,51:59] <- summary6$par
  m[i,60] <- summary6$objective*n
  
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
  if (which.min(m[i,c(10,20,30,40,50,60)])==1) {
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

# How often does generating CDM win?

n <- 1e6

l <- 1e3
m <- matrix(0,l,249)
iters1 <- rep(0,l)
iters2 <- rep(0,l)
iters3 <- rep(0,l)
iters4 <- rep(0,l)
iters5 <- rep(0,l)
iters6 <- rep(0,l)
iters7 <- rep(0,l)
iters8 <- rep(0,l)
iters9 <- rep(0,l)
iters10 <- rep(0,l)
iters11 <- rep(0,l)
iters12 <- rep(0,l)
iters13 <- rep(0,l)
iters14 <- rep(0,l)
iters15 <- rep(0,l)
iters16 <- rep(0,l)
iters17 <- rep(0,l)
iters18 <- rep(0,l)
iters19 <- rep(0,l)
iters20 <- rep(0,l)
iters21 <- rep(0,l)
iters22 <- rep(0,l)
iters23 <- rep(0,l)
iters24 <- rep(0,l)
iters25 <- rep(0,l)
iters26 <- rep(0,l)
iters27 <- rep(0,l)

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
    summary1 <- nlminb(start=start11,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='()')
  } else {
    summary1 <- nlminb(start=start21,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='()')
  }
  
  if (like12<=like22) {
    summary2 <- nlminb(start=start12,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(23)')
  } else {
    summary2 <- nlminb(start=start22,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(23)')
  }
  
  if (like13<=like23) {
    summary3 <- nlminb(start=start13,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(234)')
  } else {
    summary3 <- nlminb(start=start23,objective=CDM1logLRcpp,gradient=CDM1grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,5)),upper=rep(1,6),data=data,permutation='(234)')
  }
  
  summary4 <- nlminb(start=CDM2startRcpp(data,'()',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='()')
  summary5 <- nlminb(start=CDM2startRcpp(data,'(34)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(34)')
  summary6 <- nlminb(start=CDM2startRcpp(data,'(23)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(23)')
  summary7 <- nlminb(start=CDM2startRcpp(data,'(243)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(243)')
  summary8 <- nlminb(start=CDM2startRcpp(data,'(234)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(234)')
  summary9 <- nlminb(start=CDM2startRcpp(data,'(24)',1e-2),objective=CDM2logLRcpp,gradient=CDM2grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,6)),upper=rep(1,7),data=data,permutation='(24)')
  
  summary10 <- nlminb(start=CDM3startRcpp(data,'()',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='()')
  summary11 <- nlminb(start=CDM3startRcpp(data,'(34)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(34)')
  summary12 <- nlminb(start=CDM3startRcpp(data,'(23)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(23)')
  summary13 <- nlminb(start=CDM3startRcpp(data,'(243)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(243)')
  summary14 <- nlminb(start=CDM3startRcpp(data,'(234)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(234)')
  summary15 <- nlminb(start=CDM3startRcpp(data,'(24)',1e-2),objective=CDM3logLRcpp,gradient=CDM3grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,7)),upper=rep(1,8),data=data,permutation='(24)')
  
  summary16 <- nlminb(start=CDM4startRcpp(data,'()',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='()')
  summary17 <- nlminb(start=CDM4startRcpp(data,'(34)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(34)')
  summary18 <- nlminb(start=CDM4startRcpp(data,'(23)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(23)')
  summary19 <- nlminb(start=CDM4startRcpp(data,'(243)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(243)')
  summary20 <- nlminb(start=CDM4startRcpp(data,'(234)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(234)')
  summary21 <- nlminb(start=CDM4startRcpp(data,'(24)',1e-2),objective=CDM4logLRcpp,gradient=CDM4grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,8)),upper=rep(1,9),data=data,permutation='(24)')
  
  summary22 <- nlminb(start=CDM5startRcpp(data,'()',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='()')
  summary23 <- nlminb(start=CDM5startRcpp(data,'(34)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(34)')
  summary24 <- nlminb(start=CDM5startRcpp(data,'(23)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(23)')
  summary25 <- nlminb(start=CDM5startRcpp(data,'(243)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(243)')
  summary26 <- nlminb(start=CDM5startRcpp(data,'(234)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(234)')
  summary27 <- nlminb(start=CDM5startRcpp(data,'(24)',1e-2),objective=CDM5logLRcpp,gradient=CDM5grRcpp,control=list(eval.max=1e2,iter.max=1e2,rel.tol=1e-3/n),lower=c(-1,rep(0,9)),upper=rep(1,10),data=data,permutation='(24)')
  
  m[i,1:6] <- summary1$par
  m[i,7] <- summary1$objective*n+6*log(n)
  
  m[i,8:13] <- summary2$par
  m[i,14] <- summary2$objective*n+6*log(n)
  
  m[i,15:20] <- summary3$par
  m[i,21] <- summary3$objective*n+6*log(n)
  
  m[i,22:28] <- summary4$par
  m[i,29] <- summary4$objective*n+7*log(n)
  
  m[i,30:36] <- summary5$par
  m[i,37] <- summary5$objective*n+7*log(n)
  
  m[i,38:44] <- summary6$par
  m[i,45] <- summary6$objective*n+7*log(n)
  
  m[i,46:52] <- summary7$par
  m[i,53] <- summary7$objective*n+7*log(n)
  
  m[i,54:60] <- summary8$par
  m[i,61] <- summary8$objective*n+7*log(n)
  
  m[i,62:68] <- summary9$par
  m[i,69] <- summary9$objective*n+7*log(n)
  
  m[i,70:77] <- summary10$par
  m[i,78] <- summary10$objective*n+8*log(n)
  
  m[i,79:86] <- summary11$par
  m[i,87] <- summary11$objective*n+8*log(n)
  
  m[i,88:95] <- summary12$par
  m[i,96] <- summary12$objective*n+8*log(n)
  
  m[i,97:104] <- summary13$par
  m[i,105] <- summary13$objective*n+8*log(n)
  
  m[i,106:113] <- summary14$par
  m[i,114] <- summary14$objective*n+8*log(n)
  
  m[i,115:122] <- summary15$par
  m[i,123] <- summary15$objective*n+8*log(n)
  
  m[i,124:132] <- summary16$par
  m[i,133] <- summary16$objective*n+9*log(n)
  
  m[i,134:142] <- summary17$par
  m[i,143] <- summary17$objective*n+9*log(n)
  
  m[i,144:152] <- summary18$par
  m[i,153] <- summary18$objective*n+9*log(n)
  
  m[i,154:162] <- summary19$par
  m[i,163] <- summary19$objective*n+9*log(n)
  
  m[i,164:172] <- summary20$par
  m[i,173] <- summary20$objective*n+9*log(n)
  
  m[i,174:182] <- summary21$par
  m[i,183] <- summary21$objective*n+9*log(n)
  
  m[i,184:193] <- summary22$par
  m[i,194] <- summary22$objective*n+10*log(n)
  
  m[i,195:204] <- summary23$par
  m[i,205] <- summary23$objective*n+10*log(n)
  
  m[i,206:215] <- summary24$par
  m[i,216] <- summary24$objective*n+10*log(n)
  
  m[i,217:226] <- summary25$par
  m[i,227] <- summary25$objective*n+10*log(n)
  
  m[i,228:237] <- summary26$par
  m[i,238] <- summary26$objective*n+10*log(n)
  
  m[i,239:248] <- summary27$par
  m[i,249] <- summary27$objective*n+10*log(n)
  
  iters1[i] <- summary1$iterations
  iters2[i] <- summary2$iterations
  iters3[i] <- summary3$iterations
  iters4[i] <- summary4$iterations
  iters5[i] <- summary5$iterations
  iters6[i] <- summary6$iterations
  iters7[i] <- summary7$iterations
  iters8[i] <- summary8$iterations
  iters9[i] <- summary9$iterations
  iters10[i] <- summary10$iterations
  iters11[i] <- summary11$iterations
  iters12[i] <- summary12$iterations
  iters13[i] <- summary13$iterations
  iters14[i] <- summary14$iterations
  iters15[i] <- summary15$iterations
  iters16[i] <- summary16$iterations
  iters17[i] <- summary17$iterations
  iters18[i] <- summary18$iterations
  iters19[i] <- summary19$iterations
  iters20[i] <- summary20$iterations
  iters21[i] <- summary21$iterations
  iters22[i] <- summary22$iterations
  iters23[i] <- summary23$iterations
  iters24[i] <- summary24$iterations
  iters25[i] <- summary25$iterations
  iters26[i] <- summary26$iterations
  iters27[i] <- summary27$iterations
}

end.time <- Sys.time()
time.taken <- end.time-start.time

iters <- matrix(0,27,3)
iters[1,] <- c(min(iters1),mean(iters1),max(iters1))
iters[2,] <- c(min(iters2),mean(iters2),max(iters2))
iters[3,] <- c(min(iters3),mean(iters3),max(iters3))
iters[4,] <- c(min(iters4),mean(iters4),max(iters4))
iters[5,] <- c(min(iters5),mean(iters5),max(iters5))
iters[6,] <- c(min(iters6),mean(iters6),max(iters6))
iters[7,] <- c(min(iters7),mean(iters7),max(iters7))
iters[8,] <- c(min(iters8),mean(iters8),max(iters8))
iters[9,] <- c(min(iters9),mean(iters9),max(iters9))
iters[10,] <- c(min(iters10),mean(iters10),max(iters10))
iters[11,] <- c(min(iters11),mean(iters11),max(iters11))
iters[12,] <- c(min(iters12),mean(iters12),max(iters12))
iters[13,] <- c(min(iters13),mean(iters13),max(iters13))
iters[14,] <- c(min(iters14),mean(iters14),max(iters14))
iters[15,] <- c(min(iters15),mean(iters15),max(iters15))
iters[16,] <- c(min(iters16),mean(iters16),max(iters16))
iters[17,] <- c(min(iters17),mean(iters17),max(iters17))
iters[18,] <- c(min(iters18),mean(iters18),max(iters18))
iters[19,] <- c(min(iters19),mean(iters19),max(iters19))
iters[20,] <- c(min(iters20),mean(iters20),max(iters20))
iters[21,] <- c(min(iters21),mean(iters21),max(iters21))
iters[22,] <- c(min(iters22),mean(iters22),max(iters22))
iters[23,] <- c(min(iters23),mean(iters23),max(iters23))
iters[24,] <- c(min(iters24),mean(iters24),max(iters24))
iters[25,] <- c(min(iters25),mean(iters25),max(iters25))
iters[26,] <- c(min(iters26),mean(iters26),max(iters26))
iters[27,] <- c(min(iters27),mean(iters27),max(iters27))

count <- 0

for (i in 1:l) {
  if (which.min(m[i,c(7,14,21,29,37,45,53,61,69,78,87,96,105,114,123,133,143,153,163,173,183,194,205,216,227,238,249)])%in%c(1,4,5,10,11,16,17,22,23)) {
    count <- count+1
  }
}

iters
time.taken
count/l

# Correct CDM

# CDM 1: Wins 97.8% of time, 6.31s
# CDM 2: Wins 95.6% of time, 6.71s
# CDM 3: Wins 43.5% of time, CDM 2 or CDM 3 with correct leaf labeling wins 87.9% of time, 6.48s
# CDM 4: Wins 83.8% of time, CDM 4 or CDM 5 with correct leaf labeling wins 86.8% of time, 6.64s
# CDM 5: Wins 34.2% of time, CDM 4 or CDM 5 with correct leaf labeling wins 80.5% of time, 6.83s

# Correct topology

# CDM 1: 99.6%
# CDM 2: 99.6%
# CDM 3: 98.4%
# CDM 4: 99.6%
# CDM 5: 98.6%
