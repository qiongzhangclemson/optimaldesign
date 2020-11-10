#########################################
#Code for the experiments in Section 4.2
#########################################

rm(list=ls())

library("tidyverse")
library(gurobi)
library(xtable)
source("opt.design.r")
source("evaluation.r")
set.seed(10)

ns <- c(30) # number of patients
ps <- c(5) # number of covariates

alpha <- 0.05
betas <- matrix(0, max(ps), max(ps))
for(p in 1:max(ps)) {
  betas[1:p,p] <- 1
}
sigratios <- col_sums(betas)

betas <- rbind(1, betas)


Hall <- matrix(sample(c(-1, 1), max(ns)*(max(ps)-1), replace=TRUE), max(ns), max(ps)-1) 
Hall <- cbind(1, Hall)

Zsall <- matrix(sample(c(-1, 1), 1000*(max(ps)-1), replace=TRUE), 
             1000, max(ps)-1) 
Zsall <- cbind(1, Zsall)

re.all <- NULL
re.all.prob <- NULL
re.all.power <- NULL

for(i in 1:length(ns)) {
  for(j in 1:length(ps)) {
    n <- ns[i]
    p <- ps[j]
    
    H <- Hall[1:n, 1:p]
    Zs <- Zsall[,1:p] 

    x <- rep(c(-1, 1), each=n/2)
    nrep <- 1000
    X <- replicate(nrep, sample(x))
    vars.mean <- est(H, X, Zs)
    
    # optimal LB
    X <- opt.design.approx(H, balance=TRUE)
    vars.opt.approx <- est(H, X, Zs)
    
    # bhat
    X <- opt.design.approx1(H, balance=TRUE)
    vars.opt.add <- est(H, X, Zs)
    
    # exact
    X <- opt.design.iter(H, balance=TRUE)
    vars.opt <- est(H, X, Zs)
    
    beta0 <- betas[1:p,]
    true_sigs <- Zs %*% beta0
    
    
    re.current <- cbind(vars.mean, vars.opt.add, vars.opt.approx, vars.opt) 
    
    re.prob <- NULL
    re.power <- NULL
    for(ll in 1:length(sigratios)) {
      dev <- true_sigs[,ll]/sqrt(re.current)
      zalpha <- qnorm(1-alpha)
      re.prob0 <- cbind(sigratios[ll], pnorm(abs(dev)))
      re.power0 <- cbind(sigratios[ll],               
                        1-pnorm(zalpha-dev)+pnorm(-zalpha-dev))
      re.prob <- rbind(re.prob, re.prob0)
      re.power <- rbind(re.power, re.power0)
    }
    
    re.all <- rbind(re.all, cbind(n,p, re.current))
    re.all.prob <- rbind(re.all.prob, cbind(n,p, re.prob))
    re.all.power <- rbind(re.all.power , cbind(n,p, re.power))
  }
}

re.all.prob <- data.frame(re.all.prob)
names(re.all.prob) <- c("n","p","ratio","mean_rand","opt_add","opt_lb","opt")

re.all.power <- data.frame(re.all.power)
names(re.all.power) <- c("n","p","ratio","mean_rand","opt_add","opt_lb","opt")

re.all <- data.frame(re.all)
names(re.all) <- c("n","p","mean_rand","opt_add","opt_lb","opt")


save(re.all, re.all.prob, re.all.power,
     file=paste("n = ", n, "ill_new.RData"))
