X_t_Vasicek <- function(pi_t, t, t_init=0, x_init, mu, sigma,r, maturity){
  set.seed(123456)
  h <- t/200
  res <- x_init
  for(i in 0:199){
    res <- res+ -r*(-(mu-r)*pi_t/r-res)*h+sigma*pi_t*sqrt(h)*rnorm(1)
  }
  return(res)
}

### formulate risk constraint
qL_t_Vasi <- function(pi_t, alpha, r, mu, sigma, delta, maturity){
  x_norm <- qnorm(alpha)
  x_Lt <- -(mu-r)*(exp(r*delta)-1)*pi_t/r-sigma*sqrt((exp(2*r*delta)-1)/(2*r))*x_norm*abs(pi_t)
  #x_Lt <- (exp(-sigma * sqrt(maturity / 2) * x_norm - 
  #               (mu - sigma^2 / 2) * maturity / 2) - exp(r * delta)) * (- pi_t)
  return(x_Lt)
}


risk_con_Vasi <- function(R, alpha, r, mu, sigma, delta, maturity){
  require(rootSolve)
  solver <- function(pi_t, R, alpha, r, mu, sigma, delta, maturity) qL_t_Vasi(pi_t, alpha, r, mu, sigma, delta, maturity) - R
  result <- uniroot.all(solver, lower = -1e6, upper = 1e6, alpha = alpha, R=R, r=r, mu=mu, sigma=sigma, delta=delta, maturity=maturity)
  res <- result
  return(res)
}


##ES
qL_t_Vasi <- function(pi_t, alpha, r, mu, sigma, delta, maturity){
  x_norm <- qnorm(alpha)
  X_Lt <- -(mu-r)*(exp(r*delta)-1)*pi_t/r+sigma*sqrt((exp(2*r*delta)-1)/(2*r))*dnorm(x_norm)/alpha*abs(pi_t)
  return(X_Lt)
}

#forward simulation for Vasi
forward_sim_Vasi <- function(checking_preiods, x_0, N, R, alpha, s0, r, mu, sigma, delta, maturity){
  set.seed(123)
  
  condition <- risk_con_Vasi(R, alpha, r, mu, sigma, delta, maturity)
  condition_sim <- matrix(runif(N*(checking_preiods+2),min=-condition[1],max=condition[2]),
                          nrow=N, ncol = checking_preiods+2)
  #pi_t, t, t_init=0, x_init, mu, sigma,r, maturity
  X_sim <- matrix(runif(N*(checking_preiods+2), 
                        min=X_t_Vasicek(condition[1],maturity,t_init=0,x_0, mu,sigma, r,maturity),
                        max=X_t_Vasicek(condition[2],maturity,t_init=0,x_0, mu,sigma,r,maturity)),
                  nrow = N, ncol = checking_preiods+2)
  X_tilde_sim <- matrix(0,nrow=N, ncol=checking_preiods+2)
  X_tilde_sim[,1] <- x_0
  
  dt <- maturity/(checking_preiods+1)
  t <- 0
  for(i in 0:checking_preiods){
    t_init <- t
    t <- t_init + dt
    #print(t_init)
    #print(t)
    passing_value <- data.frame("x"=X_tilde_sim[,i+2], "X_t"=X_sim[,i+1], "pi_t"=condition_sim[,i+1])
    #pi_t, t, t_init=0, x_init, mu, sigma,r, maturity
    #print(passing_value)
    X_tilde_sim[,i+2] <- t(apply(passing_value, MARGIN = 1, function(x) x["x"]+X_t_Vasicek(x["pi_t"],t, t_init, x["X_t"], 
                                                                                           mu, sigma, r,maturity)))
  }
  res <- list("X_sim"=X_sim, "X_tilde"=X_tilde_sim, "pi_tilde"=condition_sim)
  return(res)
}




#forward simulation
forward_sim_Vasi<- function(checking_preiods, x_0, N, R, alpha, s0, r, mu, sigma, delta, maturity){
  #set.seed(123456)
  looping_periods <- ceiling(maturity/delta)
  condition <- risk_con_Vasi(R, alpha, r, mu, sigma, delta, maturity)
  #print(condition)
  if(checking_preiods==0){
    i_s <- 0
  }else{
    n_s <- looping_periods/(checking_preiods+1)
    i_s <- ceiling(seq(1,checking_preiods,by=1)*n_s)
  }
  #runif(ceiling(N*(looping_periods)),min=-1e6,max=1e6)
  #rnorm(ceiling(N*(looping_periods)),mean = 0, sd=1e3)
  condition_sim <- matrix(rnorm(ceiling(N*(looping_periods)),mean = 0, sd=1e3),
                          nrow=N, ncol = looping_periods)
  #runif(N*length(i_s),min=-3e3,max=condition)
  #max(rnorm(N*length(i_s), mean=0, sd=1e3),condition)


  
  X_sim <- matrix(rnorm(ceiling(N*(looping_periods)),mean = 0, sd=1e3),
                  nrow = N, ncol = looping_periods)
  
  #runif(ceiling(N*length(i_s)),min=-3e3, 
  #      max=X_t(condition,maturity,t_init=0,x_0, r,log_st,maturity))
  
  #max(rnorm(ceiling(N*length(i_s)),mean = 0, sd=1e3),
  #    X_t(condition,maturity,t_init=0,x_0, r,log_st,maturity))
  if(checking_preiods!=0){
    condition_sim[,i_s] <- runif(N*length(i_s),min=condition[1],max=condition[2])
    X_sim[,i_s] <-  runif(ceiling(N*length(i_s)), 
                          min=X_t_Vasicek(condition[1],maturity,t_init=0,x_0, mu,sigma, r,maturity),
                          max=X_t_Vasicek(condition[2],maturity,t_init=0,x_0, mu,sigma,r,maturity))
  }


  X_tilde_sim <- matrix(0,nrow=N, ncol=looping_periods)
  X_tilde_sim[,1] <- x_0
  dt <- delta
  t <- 0
  require(parallel)
  cl <- makeCluster(detectCores()-1)
  setDefaultCluster(cl=cl)

  for(i in 1:(looping_periods-1)){
    t_init <- t
    t <- t_init + dt
    passing_value <- data.frame("x"=X_tilde_sim[,i+1], "X_t"=X_sim[,i], "pi_t"=condition_sim[,i])
    clusterExport(cl=cl, varlist =c("checking_preiods", "x_0", "N", "R", "alpha", 
                                    "s0", "r", "mu", "sigma", "delta", "maturity", 
                                    "passing_value","t", "t_init", "X_t_Vasicek","exp_fun" ),
                        envir=environment())
    X_tilde_sim[,i+1] <- t(parApply(cl=cl, passing_value, MARGIN = 1, 
                                    function(x) x["x"]+X_t_Vasicek(x["pi_t"],t, t_init, x["X_t"],
                                                                   mu, sigma, r,maturity)))
  }
  stopCluster(cl)
  res <- list("X_sim"=X_sim, "X_tilde"=X_tilde_sim, "pi_tilde"=condition_sim)
  return(res)
}



#Backward Solution Parallel2 _version1 for decode and backup
Backward_sol_par2_Vasi <- function(dimensions, checking_preiods, x0, N, R, alpha, s0, r, mu, sigma, delta, maturity, utility){
  forward_res <- forward_sim_Vasi(checking_preiods, x0, N, R, alpha, s0, r, mu, sigma, delta, maturity)
  X_sim <- forward_res$X_sim
  X_tilde <- forward_res$X_tilde
  pi_tilde <- forward_res$pi_tilde
  condition <- risk_con_Vasi(R, alpha, r, mu, sigma, delta, maturity)
  
  column_num <- checking_preiods+2
  
  H_inverse <- inverse(utility)
  
  dt <- maturity/(checking_preiods+1)
  t <- maturity
  for(i in 0:(checking_preiods+1)){
    t <- t-dt*i
    if(t == maturity){
      V_hat_X_tilde<- apply(t(X_tilde[,checking_preiods+2]),MARGIN = 2,utility)
      #print(cat("maturity",V_hat_X_tilde))
    }else if(t<=maturity){
      #print(cat("maturity-1",V_hat_X_tilde))
      fitted_value <- apply(t(V_hat_X_tilde),MARGIN = 2, function(x) H_inverse(x)$root)
      #print(cat("fitted_value",fitted_value))
      fitting_para_X <- X_sim[,column_num-i]
      
      #print(cat("fitting_paraX",fitting_para_X))
      
      fitting_para_pi <- pi_tilde[,column_num-i]
      
      #print(cat("fitting_pi",fitting_para_pi))
      
      estimate_para <- X_tilde[,column_num-i]
      
      fitting_model <- lm(fitted_value~poly(fitting_para_X,degree = dimensions)+
                            poly(fitting_para_pi, degree = dimensions))
      
      #print(cat("model_coef", fitting_model$coefficients))
      
      epsilon <- fitting_model$residuals
      
      
      #print(cat("V_hat_X_tilde",V_hat_X_tilde))
      
      #print(cat("epsilon",epsilon))
      
      residual_model <- lm(log(epsilon^2)~poly(fitting_para_X,degree = dimensions)+
                             poly(fitting_para_pi,degree = dimensions))
      #residual_model <- 0
      pre_calc <- epsilon/exp(predict(residual_model))
      #pre_calc <-0
      #X_t by column
      
      #exp(predict(residual_model, data.frame(fitting_para=X_tilde[,column_num-i])))
      Pi_hat <- function(pi_t, X_t, pre_calc, utility){
        middle_calc <- predict(fitting_model,data.frame(c(fitting_para_X=X_t, fitting_para_pi=pi_t)))+
          #exp(predict(residual_model, data.frame(c(fitting_para_X=X_t, fitting_para_pi=pi_t))))*
          pre_calc
        middle_vector <- sapply(middle_calc, utility)
        res <- mean(middle_vector)
        return(res)
      }
      
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
      #optimise_model <- apply(t(estimate_para), 2, 
      #                        function(x) -optimParallel(condition-2, Pi_hat, lower = -Inf, upper = condition,
      #                        X_t=x, pre_calc=pre_calc, utility=utility)$value)
      optimise_model <- parApply(cl=cl, t(estimate_para), MARGIN=2, 
                                 function(x) optim(condition[2]-1, Pi_hat, lower =condition[1], upper = condition[2],
                                                   control=c(fnscale=-1),
                                                   X_t=x, pre_calc=pre_calc, utility=utility)$value)
      
      
      stopCluster(cl)
      #optimise_model <- apply(t(estimate_para), 2, 
      #                        function(x) -optim(condition-2, Pi_hat, lower = -Inf, upper = condition,
      #                                           method = "L-BFGS-B",
      #                                           X_t=x, pre_calc=pre_calc, utility=utility)$value)
      V_hat_X_tilde <- optimise_model
    }
  }
  return(mean(V_hat_X_tilde))
}


#Backward Solution Parallel2 _version2
Backward_sol_par2_Vasi <- function(dimensions, checking_preiods, x0, N, R, alpha, s0, r, mu, sigma, delta, maturity, utility){
  forward_res <- forward_sim_Vasi(checking_preiods, x0, N, R, alpha, s0, r, mu, sigma, delta, maturity)
  X_sim <- forward_res$X_sim
  X_tilde <- forward_res$X_tilde
  #print(X_tilde[1,])
  pi_tilde <- forward_res$pi_tilde
  condition <- risk_con_Vasi(R, alpha, r, mu, sigma, delta, maturity)
  
  looping_periods <- ceiling(maturity/delta)
  column_num <- looping_periods
  if(checking_preiods==0){
    i_s <- -1
  }else{
    n_s <- looping_periods/(checking_preiods+1)
    i_s <- ceiling(seq(1,checking_preiods,by=1)*n_s)
  }

  print(cat("i_s",i_s))
  print(cat("checking_period",checking_preiods))
  
  H_inverse <- inverse(utility)
  
  dt <- delta
  t <- maturity
  for(i in 0:(looping_periods-1)){
    t <- t-dt*i
    print(cat(i, " / ", (looping_periods-1)))
    if(t == maturity){
      V_hat_X_tilde<- apply(t(X_tilde[,looping_periods]),MARGIN = 2,utility)
    }else if(t<maturity){
      fitted_value <- apply(t(V_hat_X_tilde),MARGIN = 2, function(x) H_inverse(x)$root)
      fitting_para_X <- X_sim[,column_num-i]
      fitting_para_pi <- pi_tilde[,column_num-i]
      estimate_para <- X_tilde[,column_num-i]
      
      fitting_model <- lm(fitted_value~poly(fitting_para_X,degree = dimensions)+
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
        optimise_model <- optim(condition[2]-1, Pi_hat, 
                                lower = condition[1], upper = condition[2],
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
                                        "X_tilde","X_sim","forward_res","exp_fun"),
                      envir=environment())
        if(sum((column_num-i)==i_s)>0){
          optimise_model <- parApply(cl=cl, t(estimate_para), 2, 
                                    function(x) optim(condition[2]-1, Pi_hat, lower = condition[1]+100, upper = condition[2],
                                                      method="BFGS" ,control=list(fnscale=-1),
                                                      X_t=x, pre_calc=pre_calc, utility=utility)$value)
            
          
        }else{
          optimise_model <- parApply(cl=cl, t(estimate_para), MARGIN=2, 
                                     function(x) optim(condition[2]-1, Pi_hat, lower = -100, upper = 1e3,
                                                       method="BFGS",control=list(fnscale=-1),
                                                       X_t=x, pre_calc=pre_calc, utility=utility)$value)
        }
        stopCluster(cl)
      }
      
      V_hat_X_tilde <- optimise_model
    }
  }
  return(V_hat_X_tilde)
}



library(tictoc)
rec_Vasi <- NULL
for(i in 1:20){
  tic()
  rec_Vasi[i] <- Backward_sol_par2_Vasi(2, i, x_0, 1500, R, alpha, s0, r, mu, sigma, delta, maturity, S_utility_1)
  print(rec_Vasi[i])
  toc()
  print("_________________________________________")
}
plot(rec_Vasi,type = "l")


