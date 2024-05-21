
## dose-response: linear
## estimate the parameters
## poisson
## for 10.76.1.16 server

rm(list=ls())
library(dplyr)


set.seed(11052023)
## load data file
load("/home/lee9j2/ERR/datause_11052023.rda")


log.lkh <- function(data.input, x.input, x2.input, par.input){
  
  
  py <- data.input$pyr/10000
  y <- data.input$solid
  
  # sex specific
  
  male <- data.input$male
  female <- data.input$female
  
  dose <- data.input$dgy
  mdose <- male*dose
  fdose <- female*dose
  
  
  # par.input 1~4: alpha.m
  # par.input 5~8: alpha_f
  # par.input 9: theta
  # par.input 10,11: hiroshima / nagasaki
  # par.input 12,13: delta1 / delta2
  # par.input 14: phi
  # par.input 15, 16: beta.m / beta.f
  
  
  
  par1 <- par.input[1:12]
  
  e.eta <- exp(apply(x.input * par1, 2, sum))
  
  rho.d <- par.input[17] * mdose + par.input[18] * fdose
  
  par2 <- par.input[13:16]
  
  eps.z <- exp(apply(x2.input * par2, 2, sum))
  
  
  lambda <- py * e.eta * (1 + rho.d * eps.z)
  
  if(sum(lambda < 0) >= 1) log.lkh.out <- 'fail'
  else   log.lkh.out <- sum(dpois(y, lambda, log=T))
  
  return(log.lkh.out)
  
}

# test
# log.lkh(data, x, x2, par.test)

## remove some zero


byr15 <- (floor(data$year-data$age)-1915)/10
m.byr15 <- byr15 * data$male
f.byr15 <- byr15 * data$female

x <- rbind(data$male, data$mlage70, data$mlage70sq, data$mlage70qsp,
           data$female, data$flage70, data$flage70sq, data$flage70qsp
           ,m.byr15, f.byr15, data$hnic, data$nnic
)

kerma <- data$un4gy # 1 if K>4 0 otherwise
kerma[kerma==0] <- 2
kerma[kerma==1] <- 0
kerma[kerma==2] <- 1


x2 <- rbind(data$lage70 * data$male, data$lage70 * data$female , data$ex30, kerma)


## save
iter <- 5000

par.save <- matrix(0, nrow = iter, ncol = 18)
par.save[1, ] <- par <- runif(18)

for(it in 2:iter){
  
  for(p in 1:16){
    
    p.q <- rnorm(1, mean = par[p], sd = 1)
    par.q <- par; par.q[p] <- p.q
    
    num <- log.lkh(data, x, x2, par.q ) + dnorm(p.q, 0, 100^2, log=T) + dnorm(par[p], p.q, 1)
    
    den <- log.lkh(data, x, x2, par) + dnorm(par[p], 0, 100^2, log=T) + dnorm(p.q, par[p], 1)
    
    if(log(runif(1)) < (num-den)) par[p] <- p.q
    
  }
  
  # if lambda < 0 : reject
  for(p in 17:18){
    
    p.q <- rnorm(1, mean = par[p], sd = 1)
    par.q <- par; par.q[p] <- p.q
    
    if(log.lkh(data, x, x2, par.q)!="fail"){
      
      num <- log.lkh(data, x, x2, par.q ) + dnorm(p.q, 0, 100^2, log=T) + dnorm(par[p], p.q, 1)
      
      den <- log.lkh(data, x, x2, par) + dnorm(par[p], 0, 100^2, log=T) + dnorm(p.q, par[p], 1)
      
      if(log(runif(1)) < (num-den)) par[p] <- p.q
      
    }
    
  }
  
  
  par.save[it,] <- par
  
  if(it %% 500 == 0 ) save(par.save, file = paste0("/home/lee9j2/ERR/result_1206_poi", it, ".rda"))
  
  print(it)
  
  
}

save(par.save, file = "/home/lee9j2/ERR/result_total_1206_poi.rda")

