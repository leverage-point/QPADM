# n: first dimension of matrix
# p: second dimension of matrix


library(unix)
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)

library(Rcpp)
library(RcppArmadillo)
sourceCpp("arma_func.cpp")

##### update parameter #####
update_beta <- function(old_beta, old_eta, gamma, lambda, m) {
  indicate <- old_beta + old_eta / gamma
  positive_part <-
    old_beta + old_eta / gamma - lambda / (m * gamma)
  negative_part <-
    abs(old_beta + old_eta / gamma) - lambda / (m * gamma)
  new_beta <-
    (indicate>0)*pmax(positive_part,0)-(indicate<0)*pmax(negative_part,0)
  return(new_beta)
}

update_r <- function(x, y, old_u, old_beta_b, tau, gamma) {
  positive_part1 <-
    old_u / gamma + y - x %*% old_beta_b - tau / gamma
  positive_part2 <-
    -1 * old_u / gamma - y + x %*% old_beta_b + (tau - 1) / gamma
  new_r <-
    pmax(positive_part1, 0) - pmax(positive_part2, 0)
  return(new_r)
}

update_beta_b <- function(x, y, new_r, old_u, old_eta, new_beta, gamma) {
  new_beta_b <-
    (diag(ncol(x)) - t(x) %*% (solve(diag(nrow(x)) + tcrossprod(x))) %*% x) %*%
    (t(x) %*% (y - new_r + old_u / gamma) - old_eta / gamma + new_beta)
  return(new_beta_b)
}

update_u <- function(x, y, old_u, new_beta_b, new_r, gamma) {
  new_u <- old_u + gamma * (y - x %*% new_beta_b - new_r)
  return(new_u)
}

update_eta <- function(old_eta, new_beta, new_beta_b, gamma) {
  new_eta <- old_eta + gamma * (new_beta_b - new_beta)
  return(new_eta)
}
#####

object_func <- function(beta,block_param,tau,lambda,n) {
  (lambda*sum(abs(beta)) + sum(sapply(block_param, function(block_param_item) {
    sum(block_param_item$r * (tau-as.integer(block_param_item$r<0)))
  })))/n
}

QPADM <- function(x,y,tau,m,gamma,lambda,iterlim = 500,tol = 1e-4,trace = TRUE,true_beta) {
  # initialization
  n <- nrow(x)
  p <- ncol(x)
  ind <- unlist(lapply(c(1:m),function(x) rep(x,n/m)))
  
  beta <- matrix(0,nrow = p,ncol = 1)
  block_param <- list()
  for(i in 1:m) {
    beta_b_start <- rep(0,p)
    r_start <- y[(n/m*(i-1)+1):(n/m*i)]
    u_start <- rep(0,n/m)
    eta_start <- rep(0,p)
    block_param[[i]] <- list(r=r_start,beta_b=beta_b_start,u=u_start,eta=eta_start)
  }
  
  # iteration
  cl <- makeCluster(4) # detectCores() = 8
  registerDoParallel(cl,cores = 4)
  
  start_time <- proc.time() # time the algorithm
  
  precision <- 100
  block_time <- 0
  for(iternum in 1:iterlim) {
    
    loop_start_time <- proc.time()[3] # time a loop
    
    # centralization
    old_beta_mean <- as.matrix(apply(sapply(block_param, function(block_param_item) {
      block_param_item$beta_b
    }), 1, mean))
    old_eta_mean <- as.matrix(apply(sapply(block_param, function(block_param_item) {
      block_param_item$eta
    }), 1, mean))
    
    # update beta
    new_beta <- update_beta(old_beta_mean,old_eta_mean,gamma,lambda,m)
    
    # parallelization
    new_block_param <-
      foreach(b = 1:m,
              .export = c("update_r","update_beta_b","update_u","update_eta")) %dopar% {
                
                b_start_time <- proc.time()[3] # time a parallel calculation
                
                new_r <-
                  update_r(x[which(ind==b),],y[which(ind==b)],block_param[[b]]$u,block_param[[b]]$beta_b,tau,gamma)
                new_beta_b <-
                  update_beta_b(x[which(ind==b),],y[which(ind==b)],new_r,block_param[[b]]$u,block_param[[b]]$eta,new_beta,gamma)
                new_u <-
                  update_u(x[which(ind==b),],y[which(ind==b)],block_param[[b]]$u,new_beta_b,new_r,gamma)
                new_eta <-
                  update_eta(block_param[[b]]$eta,new_beta,new_beta_b,gamma)
                
                b_end_time <- proc.time()[3] # time a parallel calculation
                
                return(list(r=new_r,beta_b=new_beta_b,u=new_u,eta=new_eta,b_time=b_end_time-b_start_time))
              }
    
    max_b_time <- max(sapply(new_block_param, function(block_param_item) {
      block_param_item$b_time
    }))
    block_time <- block_time + max_b_time
    
    # check
    if(iternum > 5)
      precision <- abs((lambda*sum(abs(beta)) + sum((y - x %*% beta) * (tau - as.integer(y-x%*%beta<0))))/n - 
                         (lambda*sum(abs(new_beta)) + sum((y - x %*% new_beta) * (tau - as.integer(y-x%*%new_beta<0))))/n)
    
    beta <- new_beta
    block_param <- new_block_param
    
    AE <- sum(abs(beta-true_beta))
    LO_AE <- sum(abs(beta-true_beta)[-1])
    
    if((precision <= tol))
      break
    
    loop_end_time <- proc.time()[3] # time a loop
    if(trace==TRUE) {
      cat("iteration = ",iternum, "  max block time = ",max_b_time, 
          " loop time = ",loop_end_time-loop_start_time," precision = ",precision,
          " AE = ",AE," LO_AE = ",LO_AE,"\n")
    }
    
  }
  
  end_time <- proc.time() # time the algorithm
  
  stopCluster(cl)
  rm(cl)
  
  return(list(param = new_beta,iternum = iternum,time = end_time-start_time,block_time = block_time,precision = precision))
  
}

HBIC <- function(x,y,beta,tau) {
  n <- nrow(x)
  p <- ncol(x)
  r <- y - x %*% beta
  log(sum(r * (tau-as.integer(r<0)))) + sum(beta!=0) * log(log(n)) / n * log(p)
}

select_lambda <- function(x,y,tau,m,gamma,true_beta,start,end,step) {
  lambda_seq <- seq(start,end,by = step)
  HBIC_seq <- sapply(lambda_seq, function(lambda) {
    beta <- QPADM(x,y,tau,m,gamma,lambda,iterlim = 300,tol = 1e-5,trace = F,true_beta = true_beta)$param
    HBIC <- HBIC(x,y,beta,tau)
    AE <- sum(abs(beta-true_beta))
    LO_AE <- sum(abs(beta-true_beta)[-1])
    cat("lambda = ",lambda,", HBIC = ",HBIC,", AE = ",AE,", LO_AE = ",LO_AE,"\n")
    return(HBIC)
  })
  plot(x = lambda_seq,y = HBIC_seq,xlab = expression(lambda),ylab = "HBIC",type = 'b')
  return(lambda_seq[which.min(HBIC_seq)])
}

select_gamma <- function(x,y,tau,m,lambda,true_beta,start,end,step) {
  gamma_seq <- seq(start,end,by = step)
  HBIC_seq <- sapply(gamma_seq, function(gamma) {
    beta <- QPADM(x,y,tau,m,gamma,lambda,iterlim = 300,tol = 1e-5,trace = F,true_beta = true_beta)$param
    HBIC <- HBIC(x,y,beta,tau)
    AE <- sum(abs(beta-true_beta))
    LO_AE <- sum(abs(beta-true_beta)[-1])
    cat("gamma = ",gamma,", HBIC = ",HBIC,", AE = ",AE,", LO_AE = ",LO_AE,"\n")
    return(HBIC)
  })
  plot(x = gamma_seq,y = HBIC_seq,xlab = expression(gamma),ylab = "HBIC",type = 'b')
  return(gamma_seq[which.min(HBIC_seq)])
}


##### simulation #####
set.seed(123)
n <- 10000
p <- 6000

rlimit_as(Inf)
X <- matrix(NA,n,p)
for (j in 1:(p)) {
  X[,j] <- rnorm(n)
}
rho <- 0.5
for (j in 2:(p)){
  X[,j] <- rho*X[,j-1]+sqrt(1-rho^2)*X[,j]
}

colnames(X) <- paste("X", 1:p, sep = "")
X[, 'X1'] <- pnorm(X[, 'X1'])
epsilon <- rnorm(n, mean = 0, sd = 1)
Y <- X[, 'X6'] + X[, 'X12'] + X[, 'X15'] + X[, 'X20'] + 0.7 * epsilon * X[, 'X1']

m <- 50 # 4/10/20/50/100
tau <- 0.3 # tau = 0.3/0.5/0.7
gamma <- 0.5

true_beta <- matrix(0,nrow = p,ncol = 1)
true_beta[1] <- 0.7*qnorm(tau)
true_beta[c(6,12,15,20)] <- 1
lambda <- select_lambda(X,Y,tau,m,gamma,true_beta,start = 0.5,end = 10,step = 0.5)
gamma <- select_gamma(X,Y,tau,m,lambda,true_beta,start = 0.1,end = 1,step = 0.1)


rlimit_as(4*2^30) # bytes
rlimit_nproc(Inf) # the maximum number of processes
test <- QPADM(X,Y,tau,m,gamma=0.5,lambda=1,iterlim = 500,tol = 1e-5,true_beta = true_beta)


AE <- sum(abs(test$param-true_beta))
LO_AE <- sum(abs(test$param-true_beta)[-1])
size <- sum(test$param!=0)
P1 <- (sum(test$param[c(6,12,15,20)]!=0)==4)
P2 <- (test$param[1]!=0)
c(AE,LO_AE,size,P1,P2)















