#############################################
#The inner code called by run.experiments.R
#############################################

re.all <- NULL


for(it in 1: mirco.rep) {
    # H is covariate matrix
    H <- matrix(sample(c(-1, 1), n*(p-1), replace=TRUE), n, p-1) 
    #H <- matrix(rnorm(n*(p-1)), n, p-1) 
    H <- cbind(1, H)

    # H is covariate matrix
    Zs <- matrix(sample(c(-1, 1), 1000*(p-1), replace=TRUE), 1000, p-1) 
    #Zs <- matrix(rnorm(1000*(p-1)), 1000, p-1) 
    Zs <- cbind(1, Zs)

    
    current.time <- proc.time()
    X.opt.approx <- opt.design.approx(H, balance=TRUE)
    re <- est_lp(H, X.opt.approx, Zs)
    if(nrow(re)>0) {
      time.duration <- proc.time()-current.time
      re$method <- "LB_APPROX"
      re$IT <- it
      re$time <- time.duration[3]
    }
    re.all <- rbind(re.all, re)
    
    current.time <- proc.time()
    X.opt.approx <- opt.design.approx1(H, balance=TRUE)
    re <- est_lp(H, X.opt.approx, Zs)
    if(nrow(re)>0) {
      time.duration <- proc.time()-current.time
      re$method <- "ADDITIVE"
      re$IT <- it
      re$time <- time.duration[3]
    }
    re.all <- rbind(re.all, re)
    
    current.time <- proc.time()
    X.opt.iter <- opt.design.iter(H, balance=TRUE)
    re <- est_lp(H, X.opt.iter, Zs)
    if(nrow(re)>0) {
      time.duration <- proc.time()-current.time
      re$method <- "EXACT"
      re$IT <- it
      re$time <- time.duration[3]
    }
    re.all <- rbind(re.all, re)

    nrep <- 100

    # random balanced design
    current.time <- proc.time()
    x <- rep(c(-1, 1), each=n/2)
    X.bl <- replicate(nrep, sample(x))
    re <- est_lp(H, X.bl, Zs)
    time.duration <- proc.time()-current.time
    re$method <- "RAND"
    re$IT <- it
    re$time <- time.duration[3]
    re.all <- rbind(re.all, re)
}





re.reports <- NULL
for(it in 1:mirco.rep) {
  re.it <- re.all[re.all$IT==it,]
  re.it$time <- NULL
  random <- re.it[re.it$method=="RAND", 1]
  random.approx <- re.it[re.it$method=="RAND", 2]
  
  re.report <- re.it %>%
    filter(method !="RAND") %>%
    group_by(method, IT) %>%
    mutate(
      percentile = mean(random <= obj),
      percentile.approx = mean(random.approx <= obj.approx)
    ) %>%
    ungroup()
  re.report <- rbind(re.report, data.frame(
    obj = quantile(random, prob=c(0.01, 0.05, 0.5)),
    obj.approx = quantile(random.approx, prob=c(0.01, 0.05, 0.5)), 
    percentile = c(0.01, 0.05, 0.5),
    percentile.approx = c(0.01, 0.05, 0.5),
    method = c("RAND(1%)", "RAND(5%)", "RAND(50%)"),
    IT =it
  ))
  re.reports <- rbind(re.reports, re.report)
}
re.reports <- re.reports %>%
  mutate(
    method = factor(method, levels=c("LB_APPROX", "EXACT", "ADDITIVE", 
                                     "RAND(1%)", "RAND(5%)", "RAND(50%)")),
    IT = factor(IT)
  )

#require(lme4)
#summary(lmer(obj~method+(1|IT), re.reports))


#re.reports %>%
#  ggplot(aes(x=method, group=IT, y=obj, color=IT)) +
#  geom_line()+
#  geom_point()

re.all$n <- n
re.all$p <- p
re.reports$n <- n
re.reports$p <- p
save(re.reports, re.all, file=paste("n = ", n,"p = ",p, "re.RData"))
