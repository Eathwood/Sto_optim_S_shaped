S_utility_1 <- function(x){
  if(x>=0) return(x^0.5)
  else return(-2*abs(x)^0.5)
}

S_utility_2 <- function(x){
  if(x>=0) return(1-exp(-0.55*x))
  else return(2*(exp(0.55*x)-1))
}

inverse = function (f, lower = -1e6, upper = 1e6) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper, extendInt = "yes")[1]
}

square_inverse = inverse(function (x) x^2, 0.1, 100)
S_u1_inverse = inverse(S_utility_1)
S_u2_inverse = inverse(S_utility_2)


### Simulate log(S_t)
log_St <- function(t, s0, mu, sigma, maturity){
  res <- numeric(3)
  res[1] <- log(s0)
  drift = (mu - sigma^2) * (maturity / 2)
  res[2] <- res[1] + as.numeric(t >= maturity / 2) * (drift + sigma * sqrt(maturity/2) * rnorm(1))
  res[3] <- res[2] + as.numeric(t==maturity) * (drift + sigma * sqrt(maturity/2) * rnorm(1))
  return(res)
}

### Simulate X_t
X_t <- function(pi_t, t, t_init=0, x_init, r, log_st, maturity){ 
  dt <- t - t_init
  X_tc <- exp(r * dt) * x_init + (1 - exp(r * dt)) * pi_t
  d_log_s <- diff(log_st, 1)
  d_log_S <- sum(d_log_s[1] * (t >= maturity/2), d_log_s[2] * (t==maturity))
  res <- X_tc + pi_t * (1 - exp(-d_log_S))
  return(res)
}



### Simulate X_t with M jumps______in progress
log_St <- function(t, M, s0, mu, sigma, maturity){
  res <- numeric(M+2)
  res[1] <- log(s0)
  drift = (mu - sigma^2) * maturity/(N+1)
  for(i in 2:(M+2)){
    res[i] <- res[i-1] + as.numeric(t >= maturity/(M+1))* (drift + sigma * sqrt(maturity/(M+1)) * rnorm(1))
  }
  return(res)
}

X_t <- function(pi_t, M, t, t_init=0, x_init, r, log_st, maturity){ 
  dt <- t - t_init
  X_tc <- exp(r * dt) * x_init + (1 - exp(r * dt)) * pi_t
  d_log_s <- diff(log_st, 1)
  delta_t <- maturity/(M+1)
  d_log_S <- 0
  for(i in 1:(M+1)){
    if(t >= (delta_t*i)) d_log_S <- d_log_S+d_log_s[i]
  }
  res <- X_tc + pi_t * (1 - exp(-d_log_S))
  return(res)
}


### formulate risk constraint
qL_t <- function(pi_t, alpha, r, mu, sigma, delta, maturity){
  x_norm <- qnorm(alpha)
  x_Lt <- (exp(-sigma * sqrt(maturity / 2) * x_norm - 
                 (mu - sigma^2 / 2) * maturity / 2) - exp(r * delta)) * (- pi_t)
  return(x_Lt)
}



risk_con <- function(R, alpha, r, mu, sigma, delta, maturity){
  require(rootSolve)
  solver <- function(pi_t, R, alpha, r, mu, sigma, delta, maturity) return(qL_t(pi_t, alpha, r, mu, sigma, delta, maturity) - R)
  result <- uniroot.all(solver, lower = -1e6, upper = 1e6, alpha = alpha, R=R, r=r, mu=mu, sigma=sigma, delta=delta, maturity=maturity)
  res <- result
  return(res)
}



#forward simulation
forward_sim<- function(checking_preiods, x_0, N, R, alpha, s0, r, mu, sigma, delta, maturity,M){
  set.seed(123456)
  looping_periods <- ceiling(maturity/delta)
  condition <- risk_con(R, alpha, r, mu, sigma, delta, maturity)
  condition_normal <- R*exp(-r*delta)
  
  n_s <- looping_periods/(checking_preiods+1)
  i_s <- ceiling(seq(1,checking_preiods,by=1)*n_s)
  #runif(ceiling(N*(looping_periods)),min=-1e6,max=1e6)
  #rnorm(ceiling(N*(looping_periods)),mean = 0, sd=1e3)
  condition_sim <- matrix(rnorm(ceiling(N*(looping_periods)),mean = 0, sd=1e3),
                          nrow=N, ncol = looping_periods)
  #runif(N*length(i_s),min=-3e3,max=condition)
  #max(rnorm(N*length(i_s), mean=0, sd=1e3),condition)
  condition_sim[,i_s] <- runif(N*length(i_s),min=-3e3,max=condition_normal)
  
  X_sim <- matrix(rnorm(ceiling(N*(looping_periods)),mean = 0, sd=1e3),
                  nrow = N, ncol = looping_periods)
  
  #runif(ceiling(N*length(i_s)),min=-3e3, 
  #      max=X_t(condition,maturity,t_init=0,x_0, r,log_st,maturity))
  
  #max(rnorm(ceiling(N*length(i_s)),mean = 0, sd=1e3),
  #    X_t(condition,maturity,t_init=0,x_0, r,log_st,maturity))
  X_sim[,i_s] <-  runif(ceiling(N*length(i_s)),min=-3e3, 
                         max=X_t(condition_normal,M, maturity,t_init=0,x_0, r,
                                 log_St(maturity,M, s0, mu, sigma, maturity),
                                 maturity))
  
  if(sum(looping_periods/2==i_s)>0){
    condition_sim[,looping_periods/2] <- runif(N,min=-3e3,max=condition)
    X_sim[,looping_periods/2] <-  runif(N,min=-3e3, 
                          max=X_t(condition,M, maturity,t_init=0,x_0, r,
                                  log_St(maturity,M, s0, mu, sigma, maturity),
                                  maturity))
  }
  
  X_tilde_sim <- matrix(0,nrow=N, ncol=looping_periods)
  X_tilde_sim[,1] <- x_0
  dt <- delta
  t <- 0
  for(i in 1:(looping_periods-1)){
    t_init <- t
    t <- t_init + dt

    passing_value <- data.frame("x"=X_tilde_sim[,i+1], "X_t"=X_sim[,i], "pi_t"=condition_sim[,i])
    X_tilde_sim[,i+1] <- t(apply(passing_value, MARGIN = 1, function(x) x["x"]+X_t(x["pi_t"],M, t, t_init, x["X_t"], r,
                                                                                   log_St(maturity,M, s0, mu, sigma, maturity),
                                                                                   maturity)))

  }
  res <- list("X_sim"=X_sim, "X_tilde"=X_tilde_sim, "pi_tilde"=condition_sim)
  return(res)
}


#Backward Solution Parallel2
Backward_sol_par2<- function(dimensions, checking_preiods, x0, N, R, alpha, s0, r, mu, sigma, delta, maturity, utility,M){
  forward_res <- forward_sim(checking_preiods, x0, N, R, alpha, s0, r, mu, sigma, delta, maturity,M)
  X_sim <- forward_res$X_sim
  X_tilde <- forward_res$X_tilde
  #print(X_tilde[1,])
  pi_tilde <- forward_res$pi_tilde
  condition <- risk_con(R, alpha, r, mu, sigma, delta, maturity)

  looping_periods <- ceiling(maturity/delta)
  column_num <- looping_periods
  n_s <- looping_periods/(checking_preiods+1)
  i_s <- ceiling(seq(1,checking_preiods,by=1)*n_s)
  
  print(cat("i_s",i_s))
  print(cat("checking_period",checking_preiods))
  
  H_inverse <- inverse(utility)
  
  dt <- delta
  t <- maturity
  for(i in 0:(looping_periods-1)){
    t <- t-dt*i

    if(t == maturity){
      V_hat_X_tilde<- apply(t(X_tilde[,looping_periods]),MARGIN = 2,utility)
    }else if(t<maturity){
      fitted_value <- apply(t(V_hat_X_tilde),MARGIN = 2, function(x) H_inverse(x)$root)
      fitting_para_X <- X_sim[,column_num-i]
      fitting_para_pi <- pi_tilde[,column_num-i]
      estimate_para <- X_tilde[,column_num-i]
      

      
      fitting_model <- lm((fitted_value+1e-19)~poly(fitting_para_X,degree = dimensions)+
                            poly(fitting_para_pi, degree = dimensions))
      
      
      epsilon <- fitting_model$residuals

      residual_model <- lm(log(epsilon^2+1e-19)~poly(fitting_para_X,degree = dimensions)+
                             poly(fitting_para_pi,degree = dimensions))
      
      pre_calc <- epsilon/exp(predict(residual_model))

      
      Pi_hat <- function(pi_t, X_t, pre_calc, utility){
        middle_calc <- predict(fitting_model,data.frame(fitting_para_X=X_t, fitting_para_pi=pi_t))+
          exp(predict(residual_model, data.frame(fitting_para_X=X_t, fitting_para_pi=pi_t)))*
          pre_calc
        middle_vector <- sapply(middle_calc, utility)
        res <- mean(middle_vector[which(!is.na(middle_vector))])
        return(res)
      }
      
      if(i==looping_periods-1){
        optimise_model <- optim(condition-1, Pi_hat, 
                                #method="BFGS" ,
                                lower = -1e3, upper = 1e3,
                                control=c(fnscale=-1),
                                X_t=estimate_para[1], pre_calc=pre_calc, utility=utility)$value
      }else{
        require(parallel)
        cl <- makeCluster(detectCores()-1)
        setDefaultCluster(cl=cl)
        clusterExport(cl=cl, varlist =c("dimensions", "checking_preiods", "x0", "N", "R", "alpha", 
                                        "s0", "r", "mu", "sigma", "delta", "maturity", "utility",
                                        "Pi_hat", "pre_calc","i","residual_model","epsilon","fitting_model",
                                        "estimate_para","fitting_para_pi","fitting_para_X","fitted_value",
                                        "V_hat_X_tilde","t","dt","H_inverse","column_num","condition","pi_tilde",
                                        "X_tilde","X_sim","forward_res"),
                      envir=environment())
        if(sum((column_num-i)==i_s)>0){
          if((column_num-i)==(looping_periods/2)){
            
            optimise_model <- parApply(cl=cl,t(estimate_para), 2, 
                                    function(x) optim(condition-1, Pi_hat, 
                                                      #method="BFGS" ,
                                                      lower = -1e3, upper = condition,
                                                      control=list(fnscale=-1),
                                                      X_t=x, pre_calc=pre_calc, utility=utility)$value)
            #stopCluster(cl)
            
          }else{
            optimise_model <- parApply(cl=cl, t(estimate_para), 2, 
                                    function(x) optim(R*exp(-r*delta)-1, Pi_hat, 
                                                      #method="BFGS" ,
                                                      lower = -1e3, upper = R*exp(-r*delta),
                                                      control=list(fnscale=-1),
                                                      X_t=x, pre_calc=pre_calc, utility=utility)$value)
            #stopCluster(cl)
          }
        }else{
          optimise_model <- parApply(cl=cl, t(estimate_para), MARGIN=2, 
                                     function(x) optim(condition-1, Pi_hat, 
                                                       #method="BFGS" ,
                                                       lower = -1e3, upper = 1e3,
                                                       control=list(fnscale=-1),
                                                       X_t=x, pre_calc=pre_calc, utility=utility)$value)
          #stopCluster(cl)
        }
        stopCluster(cl)
      }
      
      V_hat_X_tilde <- optimise_model
      #print(cat(i," / ",(looping_periods-1)))
    }
  }
  return(V_hat_X_tilde)
}



library(tictoc)
alpha <- 0.01
rec_U11_001M3 <- NULL
for(i in 1:10){
  tic()
  rec_U11_001M3[i] <- Backward_sol_par2(2, i, x_0, 50, R, alpha, s0, r, mu, sigma, delta, maturity*3, S_utility_1,2)
  print(rec_U11_001M3[i])
  toc()
  print("_________________________________________")
}

