##############################################
#The experiments and results for Section 4.1
#############################################

rm(list=ls())

library("tidyverse")
library(gurobi)
library(xtable)
source("opt.design.r")
source("evaluation.r")
set.seed(10)

mirco.rep <- 5# number of mirco.replication

ns <- c(60, 100, 120, 150) # number of patients
ps <- c(4, 6, 10, 15, 20) # number of covariates



for(i in 1:length(ns)) {
  for(j in 1:length(ps)) {
    n <- ns[i]
    p <- ps[j]
    source("main.r")
    print(c(n,p))
    #save(re.all, file=paste("n = ", n,"p = ",p, "re.RData"))
  }
}

#################################################
# The following code generate figure 1, figure 5
#################################################

# present results
collect.re <- NULL
for(i in 1:length(ns)) {
  for(j in 1:length(ps)) {
    n <- ns[i]
    p <- ps[j]
    load(paste("n = ", n,"p = ",p, "re.RData"))
    collect.re <- rbind(collect.re, re.reports)
  }
}



ll <- expand.grid(paste("n =", ns), paste("p =", ps))
ll <- paste(ll$Var1, ll$Var2, sep=", ")

collect.re %>%
  mutate(
    p = paste("p =",p), 
    p = factor(p, levels=paste("p =", ps)),
    n = paste("n =",n), 
    n = factor(n, levels=paste("n =", ns)),
    np = paste(n, p, sep=", "),
    np = factor(np, levels=ll)
  ) %>% 
  ggplot(aes(x=obj, y=obj.approx)) +
  geom_point()+
  facet_wrap(~np, scale="free", ncol=4)+
  theme_bw()+
  ylab("Surrogate Objective Value")+
  xlab("Original Objective Value")

ggsave("compare1.pdf", width=8, height=10)

  


collect.re %>%
  mutate(
    p = paste("p =",p), 
    p = factor(p, levels=paste("p =", ps)),
    n = paste("n =",n), 
    n = factor(n, levels=paste("n =", ns)),
    np = paste(n, p, sep=", "),
    np = factor(np, levels=ll)
  ) %>%
  ggplot(aes(x=method, y=obj, color=IT, group=IT))+
  geom_point(aes(shape=IT))+
  geom_line(aes(linetype=IT))+
  facet_wrap(~np, scale="free_x", ncol=4)+
  coord_flip()+
  ylab("Original Objective Value")+
  theme_bw()+
  theme(legend.position="none")
ggsave("true_obj.pdf", width=8, height=10)


collect.re %>%
  mutate(
    p = paste("p =",p), 
    p = factor(p, levels=paste("p =", ps)),
    n = paste("n =",n), 
    n = factor(n, levels=paste("n =", ns)),
    np = paste(n, p, sep=", "),
    np = factor(np, levels=ll)
  ) %>%
  ggplot(aes(x=method, y=obj.approx, color=IT, group=IT))+
  geom_point(aes(shape=IT))+
  geom_line(aes(linetype=IT))+
  facet_wrap(~np, scale="free_x", ncol=4)+
  coord_flip()+
  ylab("Surrogate Objective Value")+
  theme_bw()+
  theme(legend.position="none")
ggsave("approx_obj.pdf", width=8, height=10)








