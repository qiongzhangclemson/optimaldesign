rm(list=ls())

library("tidyverse")
library(gurobi)
library(xtable)
source("opt.design.r")
source("evaluation.r")
set.seed(10)

#-----------------------------
# Results: variance reduction
#-----------------------------


ns <- c(60, 100, 120, 150) # number of patients
#ns <- 300

results <- NULL
for(i in 1:length(ns)) {
    n <- ns[i]
    load(file=paste("n = ", n, "ill_new.RData"))
    results <- rbind(results, re.all)
}

load(file=paste("n = ", 300, "ill_new.RData"))
re.all$opt <- NA
results <- rbind(results, re.all)


names(results)[4:6] <- c("ADDITIVE","LB_APPROX", "EXACT") 
results %>%
  gather(method, var ,-n, -p, -mean_rand) %>%
  mutate(
    var_re=(mean_rand-var)/mean_rand,
    p = factor(p)
  ) %>%
  ggplot(aes(x=p, y=var_re))+
  geom_boxplot()+
  facet_grid(method~n, scale="free_x", labeller = label_both)+
  ylab("Variance Reduction")+
  theme_bw()



ggsave("var_re.pdf", width=7.5, height=4.5)




#-------------------------------------------
# Results: Probability of Correct Selection Improvement
#-------------------------------------------


ns <- c(60, 100, 150) # number of patients
#ns <- 300

results <- NULL
for(i in 1:length(ns)) {
  n <- ns[i]
  load(file=paste("n = ", n, "ill_new.RData"))
  results <- rbind(results, re.all.prob)
}


results <- results[,-4]
names(results)[4:6] <- c("ADDITIVE","LB_APPROX", "EXACT") 
results <- results %>%
  gather(method, prob ,-n, -p, -ratio, -ADDITIVE) %>%
  group_by(method, n, ratio, p) %>%
  mutate(
    prob_improv=(prob-ADDITIVE)/ADDITIVE,
    quantile = rank(-ADDITIVE)/length(ADDITIVE),
    quantile = round(quantile*5)/5
  ) 
results%>%
  filter(ratio%in%c(0.01, 0.05,0.1), p==20, method!="RAND") %>%
  ggplot(aes(x=factor(quantile), y=prob_improv, color=method))+
  geom_boxplot()+
  facet_grid(ratio~n, labeller = label_both)+
  ylab("Probability of Correct Selection Improvement")+
  xlab("Quantile of Variance")+
  theme_bw()+
  theme(legend.position = "bottom")



ggsave("prob_improve.pdf", width=8, height=8)


#-------------------------------------------
# Results: Probability of Correct Selection 
#-------------------------------------------


ns <- c(60, 100, 150) # number of patients

results <- NULL
for(i in 1:length(ns)) {
  n <- ns[i]
  load(file=paste("n = ", n, "ill_new.RData"))
  results <- rbind(results, re.all.prob)
}


results <- results[,-4]
names(results)[4:6] <- c("ADDITIVE","LB_APPROX", "EXACT") 
results <- results %>%
  gather(method, prob ,-n, -p, -ratio) %>%
  group_by(method, n, ratio, p) %>%
  mutate(
    quantile = rank(-prob)/length(prob),
    quantile = round(quantile*5)/5
  ) 
results%>%
  filter(ratio%in%c(0.01, 0.05,0.1), p==20) %>%
  ggplot(aes(x=factor(quantile), y=prob, color=method))+
  geom_boxplot()+
  facet_grid(ratio~n, labeller = label_both)+
  ylab("Probability of Correct Selection")+
  xlab("Quantile of Variance")+
  theme_bw()+
  theme(legend.position = "bottom")



ggsave("prob.pdf", width=8, height=8)