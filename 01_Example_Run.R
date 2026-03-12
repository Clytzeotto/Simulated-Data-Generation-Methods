# ================================================================================
# Script Name：01_Example_Run.R
# The Examples of Table2, Table3 and TableA15
# Package：[simsurv、progress]
# ================================================================================

library("simsurv")
library("progress")
source("D:/xxx/xxx/xxx/xxx/00_core_functions.R")

#Table2:Results of key parameters and censoring rates under different simulation conditions
#Table2:ID1
t1 <- proc.time()
thetacal_out <- thetacal(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                         gamma2=0.75,lambda2=0.5,mixp=0,
                         cens.p=0.2,size=1000,binomp=0.5,
                         beta1=0,beta2=0,maxt=6,mean=0,std=1,
                         threshold=0.00001,censdist = "E")
t2 <- proc.time()
t <- t2-t1
print(paste0('Execution time：',round(t[3][[1]],3),'s'))

theta <- thetacal_out[[1]]
show(theta)
test1 <- censoringcal(starti=2001,simusize=1000,simusizeplus=0,
                     gamma1=1.4,lambda1=0.1,
                     gamma2=0.75,lambda2=0.5,mixp=0,
                     cens.p=0.2,size=1000,binomp=0.5,
                     beta1=0,beta2=0,maxt=6,mean=0,std=1,censdist="E")
test1$theta <- theta
write.csv(test1,"D:/xxx/xxx/xxx/xxx/test1.csv")
summary(test1$CensPAll)

#Table2:ID2
t1 <- proc.time()
thetacal_out <- thetacal(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                         gamma2=0.75,lambda2=0.5,mixp=0,
                         cens.p=0.2,size=1000,binomp=0.5,
                         beta1=0.5,beta2=0,maxt=6,mean=0,std=1,threshold=0.00001,censdist = "E")
t2 <- proc.time()
t <- t2-t1
print(paste0('Execution time：',round(t[3][[1]],3),'s'))

theta <- thetacal_out[[1]]
show(theta)
test2 <- censoringcal(starti=2001,simusize=1000,simusizeplus=0,
                      gamma1=1.4,lambda1=0.1,
                      gamma2=0.75,lambda2=0.5,mixp=0,
                      cens.p=0.2,size=1000,binomp=0.5,
                      beta1=0.5,beta2=0,maxt=6,mean=0,std=1,censdist="E")
test2$theta <- theta
write.csv(test2,"D:/xxx/xxx/xxx/xxx/test2.csv")
summary(test2$CensPAll)

#Table2:ID4
t1 <- proc.time()
thetacal_out <- thetacal(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                         gamma2=0.75,lambda2=0.5,mixp=0,
                         cens.p=0.2,size=1000,binomp=0.5,
                         beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,threshold=0.00001,censdist = "E")
t2 <- proc.time()
t <- t2-t1
print(paste0('Execution time：',round(t[3][[1]],3),'s'))

theta <- thetacal_out[[1]]
show(theta)

test4 <- censoringcal(starti=2001,simusize=1000,simusizeplus=0,
                      gamma2=0.75,lambda2=0.5,mixp=0,
                      cens.p=0.2,size=1000,binomp=0.5,
                      beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,censdist="E")
test4$theta <- theta
write.csv(test4,"D:/xxx/xxx/xxx/xxx/test4.csv")
summary(test4$CensPAll)

#Table3:Comparison of parametric methods and kernel density estimation in estimating covariate effects:
#computational time, key parameters, and simulated censoring rates
#Table3:ID2
time2 <- data.frame()
temp <- c()
abT10 <- c()
pb <- progress_bar$new(
  format="[:bar] Percent of Completed Iteration: :percent   ETA: :eta  ELA: :elapsedfull",
  total=100, clear=FALSE, width=80
)
for (i in 1:100){
  t1 <- proc.time()
  thetacal_out <- thetacal(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                           gamma2=0.75,lambda2=0.5,mixp=0,
                           cens.p=0.2,size=1000,binomp=0.5,
                           beta1=0.5,beta2=0,maxt=6,mean=0,std=1,threshold=0.00001,censdist = "E")
  t2 <- proc.time()
  t <- t2-t1
  theta <- thetacal_out[[1]]
  
  t1_a <- proc.time()
  thetacal_out_a <- thetacal_a(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                               gamma2=0.75,lambda2=0.5,mixp=0,
                               cens.p=0.2,size=1000,binomp=0.5,
                               beta1=0.5,beta2=0,maxt=6,mean=0,std=1,threshold=0.00001,censdist = "E",unirootlower=0,unirootupper=50)
  t2_a <- proc.time()
  t_a <- t2_a-t1_a
  theta_a <- thetacal_out_a[[1]]
  temp <- c(2,i,thetacal_out[[1]],t[3][[1]],thetacal_out_a[[1]],t_a[3][[1]])
  
  time2 <- rbind(time2,temp)
  pb$tick()
}
colnames(time2) <- c("ID","i","theta","t","theta_a","t_a")
write.csv(time2,"D:/xxx/xxx/xxx/xxx/time2.csv")
theta <- thetacal_out_a[[1]]
show(theta)
test2_a <- censoringcal(starti=2001,simusize=1000,simusizeplus=0,
                        gamma2=0.75,lambda2=0.5,mixp=0,
                        cens.p=0.2,size=1000,binomp=0.5,
                        beta1=0.5,beta2=0,maxt=6,mean=0,std=1,censdist="E")
test2_a$theta_a <- theta
write.csv(test2_a,"D:/xxx/xxx/xxx/xxx/test2_a.csv")
summary(test2_a$CensPAll)

#Table3:ID18
time18 <- data.frame()
temp <- c()
abT10 <- c()
pb <- progress_bar$new(
  format="[:bar] Percent of Completed Iteration: :percent   ETA: :eta  ELA: :elapsedfull",
  total=100, clear=FALSE, width=80
)
for (i in 1:100){
  t1 <- proc.time()
  thetacal_out <- thetacal(seedtheta=2001,gamma1=0.8,lambda1=0.4,
                           gamma2=2,lambda2=0.1,mixp=0.75,
                           cens.p=0.2,size=1000,binomp=0.5,
                           beta1=0.5,beta2=0,maxt=6,mean=0,std=1,threshold=0.00001,censdist = "E")
  t2 <- proc.time()
  t <- t2-t1
  theta <- thetacal_out[[1]]
  
  t1_a <- proc.time()
  thetacal_out_a <- thetacal_a(seedtheta=2001,gamma1=0.8,lambda1=0.4,
                               gamma2=2,lambda2=0.1,mixp=0.75,
                               cens.p=0.2,size=1000,binomp=0.5,
                               beta1=0.5,beta2=0,maxt=6,mean=0,std=1,threshold=0.00001,censdist = "E",unirootlower=0,unirootupper=50)
  t2_a <- proc.time()
  t_a <- t2_a-t1_a
  theta_a <- thetacal_out_a[[1]]
  temp <- c(18,i,thetacal_out[[1]],t[3][[1]],thetacal_out_a[[1]],t_a[3][[1]])
  
  time18 <- rbind(time18,temp)
  pb$tick()
}
colnames(time18) <- c("ID","i","theta","t","theta_a","t_a")
write.csv(time18,"D:/xxx/xxx/xxx/xxx/time18.csv")
theta <- thetacal_out_a[[1]]
show(theta)
test18_a <- censoringcal(starti=2001,simusize=1000,simusizeplus=0,
                         gamma1=0.8,lambda1=0.4,
                         gamma2=2,lambda2=0.1,mixp=0.75,
                         cens.p=0.2,size=1000,binomp=0.5,
                         beta1=0.5,beta2=0,maxt=6,mean=0,std=1,censdist="E")
test18_a$theta_a <- theta
write.csv(test18_a,"D:/xxx/xxx/xxx/xxx/test18_a.csv")
summary(test18_a$CensPAll)

#Appendix A.TableA15:Key Parameters and Censoring Rates at Different Probability Density Thresholds
#TableA15:ID4--the threshold is 0.00001
t1 <- proc.time()
thetacal_out <- thetacal(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                         gamma2=0.75,lambda2=0.5,mixp=0,
                         cens.p=0.2,size=1000,binomp=0.5,
                         beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,threshold=0.00001,censdist = "E")
t2 <- proc.time()
t <- t2-t1
print(paste0('Execution time：',round(t[3][[1]],3),'s'))

theta <- thetacal_out[[1]]
show(theta)

test4_4 <- cc(starti=2001,simusize=1000,simusizeplus=0,
              gamma1=1.4,lambda1=0.1,
              gamma2=0.75,lambda2=0.5,mixp=0,
              cens.p=0.2,size=1000,binomp=0.5,
              beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,censdist="E")
test4_4$theta <- theta
test4_4$t <- t[3][[1]]
test4_4$id <- "4_4"
write.csv(test4_4,"D:/xxx/xxx/xxx/xxx/test4_4.csv")
summary(test4_4$CensPAll)

#TableA15:ID4--the threshold is 0.0001
t1 <- proc.time()
thetacal_out <- thetacal(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                         gamma2=0.75,lambda2=0.5,mixp=0,
                         cens.p=0.2,size=1000,binomp=0.5,
                         beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,threshold=0.0001,censdist = "E")
t2 <- proc.time()
t <- t2-t1
print(paste0('Execution time：',round(t[3][[1]],3),'s'))

theta <- thetacal_out[[1]]
show(theta)

test4_3 <- cc(starti=2001,simusize=1000,simusizeplus=0,
              gamma1=1.4,lambda1=0.1,
              gamma2=0.75,lambda2=0.5,mixp=0,
              cens.p=0.2,size=1000,binomp=0.5,
              beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,censdist="E")
test4_3$theta <- theta
test4_3$t <- t[3][[1]]
test4_3$id <- "4_3"
write.csv(test4_3,"D:/xxx/xxx/xxx/xxx/test4_3.csv")
summary(test4_3$CensPAll)

#TableA15:ID4--the threshold is 0.001
t1 <- proc.time()
thetacal_out <- thetacal(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                         gamma2=0.75,lambda2=0.5,mixp=0,
                         cens.p=0.2,size=1000,binomp=0.5,
                         beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,threshold=0.001,censdist = "E")
t2 <- proc.time()
t <- t2-t1
print(paste0('Execution time：',round(t[3][[1]],3),'s'))

theta <- thetacal_out[[1]]
show(theta)

test4_2 <- cc(starti=2001,simusize=1000,simusizeplus=0,
              gamma1=1.4,lambda1=0.1,
              gamma2=0.75,lambda2=0.5,mixp=0,
              cens.p=0.2,size=1000,binomp=0.5,
              beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,censdist="E")
test4_2$theta <- theta
test4_2$t <- t[3][[1]]
test4_2$id <- "4_2"
write.csv(test4_2,"D:/xxx/xxx/xxx/xxx/test4_2.csv")
summary(test4_2$CensPAll)

#TableA15:ID4--the threshold is 0.01
t1 <- proc.time()
thetacal_out <- thetacal(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                         gamma2=0.75,lambda2=0.5,mixp=0,
                         cens.p=0.2,size=1000,binomp=0.5,
                         beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,threshold=0.01,censdist = "E")
t2 <- proc.time()
t <- t2-t1
print(paste0('Execution time：',round(t[3][[1]],3),'s'))

theta <- thetacal_out[[1]]
show(theta)

test4_1 <- cc(starti=2001,simusize=1000,simusizeplus=0,
              gamma1=1.4,lambda1=0.1,
              gamma2=0.75,lambda2=0.5,mixp=0,
              cens.p=0.2,size=1000,binomp=0.5,
              beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,censdist="E")
test4_1$theta <- theta
test4_1$t <- t[3][[1]]
test4_1$id <- "4_1"
write.csv(test4_1,"D:/xxx/xxx/xxx/xxx/test4_1.csv")
summary(test4_1$CensPAll)

#TableA15:ID4--the threshold is 0.1
t1 <- proc.time()
thetacal_out <- thetacal(seedtheta=2001,gamma1=1.4,lambda1=0.1,
                         gamma2=0.75,lambda2=0.5,mixp=0,
                         cens.p=0.2,size=1000,binomp=0.5,
                         beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,threshold=0.1,censdist = "E")
t2 <- proc.time()
t <- t2-t1
print(paste0('Execution time：',round(t[3][[1]],3),'s'))

theta <- thetacal_out[[1]]
show(theta)

test4_0 <- cc(starti=2001,simusize=1000,simusizeplus=0,
              gamma1=1.4,lambda1=0.1,
              gamma2=0.75,lambda2=0.5,mixp=0,
              cens.p=0.2,size=1000,binomp=0.5,
              beta1=0.5,beta2=0.1,maxt=6,mean=0,std=1,censdist="E")
test4_0$theta <- theta
test4_0$t <- t[3][[1]]
test4_0$id <- "4_0"
write.csv(test4_0,"D:/xxx/xxx/xxx/xxx/test4_0.csv")
summary(test4_0$CensPAll)
