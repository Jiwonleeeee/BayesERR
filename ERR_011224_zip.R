
## dose-response: linear
## estimate the parameters - without latent variables
## update only the coefficients
## zip
## for arcc

rm(list=ls())
library(dplyr)


set.seed(11052023)
## load data file
load("/home/lee9j2/ERR/datause_11052023.rda")
# load("/Users/wonny/Documents/ERR_local/smalldata.rda")

log.lkh <- function(data.input, x.input, x2.input, par.input, B.input){
  
  
  py <- data.input$pyr/10000
  y <- data.input$solid
  
  # sex specific
  
  male <- data.input$male
  female <- data.input$female
  
  dose <- data.input$dgy
  mdose <- male*dose
  fdose <- female*dose
  
  
  par1 <- par.input[1:12]
  
  e.eta <- exp(apply(x.input * par1, 2, sum))
  
  rho.d <- par.input[17] * mdose + par.input[18] * fdose
  
  par2 <- par.input[13:16]
  
  eps.z <- exp(apply(x2.input * par2, 2, sum))
  
  lambda <- py * e.eta * (1 + rho.d * eps.z)
  
  B0 <- which(B.input==0)
  
  lambda <- lambda[B0]
  y <- y[B0]
  
  
  
  if(sum(lambda < 0) >= 1) log.lkh.out <- 'fail'
  else   log.lkh.out <- sum(dpois(y, lambda, log=T))
  
  return(log.lkh.out)
  
}


lb.ftn <- function(data.input, x.input, x2.input, par.input){
  
  
  py <- data.input$pyr/10000
  y <- data.input$solid
  n <- length(y)
  # sex specific
  
  male <- data.input$male
  female <- data.input$female
  
  dose <- data.input$dgy
  mdose <- male*dose
  fdose <- female*dose
  

  par1 <- par.input[1:12]
  
  e.eta <- exp(apply(x.input * par1, 2, sum))
  
  rho.d <- par.input[17] * mdose + par.input[18] * fdose
  
  par2 <- par.input[13:16]
  
  eps.z <- exp(apply(x2.input * par2, 2, sum))
  
  lambda <- py * e.eta * (1 + rho.d * eps.z)
  
  return(lambda)
  
}


byr15 <- (floor(data$year-data$age)-1915)/10
m.byr15 <- byr15 * data$male
f.byr15 <- byr15 * data$female

x <- rbind(data$male, data$mlage70, data$mlage70sq, data$mlage70qsp,
           data$female, data$flage70, data$flage70sq, data$flage70qsp
           ,m.byr15, f.byr15, data$hnic, data$nnic
)


z <- rbind(data$male, data$female, data$dgy, data$ex30, data$hnic, data$nnic)


kerma <- data$un4gy # 1 if K>4 0 otherwise
kerma[kerma==0] <- 2
kerma[kerma==1] <- 0
kerma[kerma==2] <- 1


x2 <- rbind(data$lage70 * data$male, data$lage70 * data$female , data$ex30, kerma)


## save
iter <- 5000

par.save <- matrix(0, nrow = iter, ncol = 24)
par.save[1, ] <- par <- runif(24)

# initial B
B <- numeric(nrow(data))
half <- sample(which(data$solid==0), length(which(data$solid==0))/2)
B[half] <- 1
B_save <- B


## initial pi
kappa <- par[19:24]

z <- rbind(data$male, data$female, data$dgy, data$ex30, data$hnic, data$nnic)

pi <- exp(apply(z* kappa, 2, sum))/(1 + exp(apply(z* kappa, 2, sum))) # initial pi

# initial lambda
lambda <- lb.ftn(data, x, x2, par)


y <- data$solid

for(it in 2:iter){
  
  
  ## gamma 
  
  for(p in 1:16){
    
    p.q <- rnorm(1, mean = par[p], sd = 1)
    par.q <- par; par.q[p] <- p.q
    
    num <- log.lkh(data, x, x2, par.q, B ) + dnorm(p.q, 0, 100^2, log=T) + dnorm(par[p], p.q, 1)
    
    den <- log.lkh(data, x, x2, par, B) + dnorm(par[p], 0, 100^2, log=T) + dnorm(p.q, par[p], 1)
    
    if(log(runif(1)) < (num-den)) par[p] <- p.q
    
  }
  
  # if lambda < 0 : reject
  for(p in 17:18){
    
    p.q <- rnorm(1, mean = par[p], sd = 1)
    par.q <- par; par.q[p] <- p.q
    
    if(log.lkh(data, x, x2, par.q, B)!="fail"){
      
      num <- log.lkh(data, x, x2, par.q, B ) + dnorm(p.q, 0, 100^2, log=T) + dnorm(par[p], p.q, 1)
      
      den <- log.lkh(data, x, x2, par, B) + dnorm(par[p], 0, 100^2, log=T) + dnorm(p.q, par[p], 1)
      
      if(log(runif(1)) < (num-den)) par[p] <- p.q
      
    }
    
  }
  
  ## B update only for y==0, 
  y0.which <- which(y==0)
  
  prob <- pi[y0.which] / (pi[y0.which] + exp(-lambda[y0.which])*(1-pi[y0.which]))
  
  B[y0.which] <- sapply(1:sum(y==0), function(x){sample(c(1,0), 1, prob=c(prob[x], (1-prob[x])))} ) 
  
  
  ## kappa: par 19~24
  for(p in 19:24){
    
    p.q <- rnorm(1, mean = par[p], sd = 1)
    par.q <- par; par.q[p] <- p.q
    
    # for each par.q - recalculate pi
    pi.q <- exp(apply(z* par.q[19:24], 2, sum))/(1 + exp(apply(z* par.q[19:24], 2, sum)))
    
    num <- sum( B * log(pi.q) + (1-B) * log(1-pi.q)   ) + dnorm(p.q, 0, 100^2, log=T) + dnorm(par[p], p.q, 1)
    
    den <- sum( B * log(pi) + (1-B) * log(1-pi)   ) + dnorm(par[p], 0, 100^2, log=T) + dnorm(p.q, par[p], 1)
    
    if(log(runif(1)) < (num-den)) par[p] <- p.q
    
    
  }
  
  
  
  
  
  
  par.save[it,] <- par
  B_save <- B_save + B
  
  if(it %% 1000 == 0 ) save(par.save, B_save, file = paste0("/home/lee9j2/ERR/result_011224_latent", it, ".rda"))
  
  print(it)
  
  
}

save(par.save,  B_save, file = "/home/lee9j2/ERR/result/result_011224_latent_total.rda")

