# Dynamic Sparse Factor Analysis: MCMC Implementation
# Call:
#   results = MCMC_DSS(Y, phi_w, phi_beta_init = 0.95, phiest = TRUE,
#     lambda0, lambda1, Theta, delta, n0, d0, RV = FALSE, NITER = 10000)

####################################################################
####################################################################

# Kalman Filter & Smoother for Latent Factors

# Input: 
#   phi_w: AR coefficient of latent factors
#   Beta: a list of T matrices; pxk matrix B_t = Beta[[t]]
#   Sigma: a list of T matrices; pxp matrix S_t = Sigma[[t]]
#   Y: pxT Data matrix

# Output:
#   result: a list with three components
#       result$w: a list of T vectors; k-vector w_{t | T} = result$w[[t]]
#       result$V: a list of T matrices; kxk matrix V_{t | T} = result$V[[t]]

# Other Notations:
#   w1: a list of T vectors: k-vector w_{t|t} = w1[[t]]
#   w12: a list of T vectors: k-vector w_{t|t-1} = w12[[t]]
#   V1: a list of T matrices: kxk matrix V_{t|t} = V1[[t]]
#   V12: a list of T matrices: kxk matrix V_{t|t-1} = V12[[t]]

#__________________________________________________________________

Kalman_FS = function(phi_w, Beta, Sigma, Y){
  
  T = length(Beta)
  p = nrow(Beta[[1]])
  k = ncol(Beta[[1]])
  
# Initialize w_0|0 and V_0|0
  w_0 = rep(0,k)
  V_0 = 1/(1-phi_w^2)*diag(k)
  
  w1 = list()
  w12 = list()
  V1 = list()
  V12 = list()
  
  result = list()
  w = list()
  V = list()
  
# Forward Prediction & Correction Step for t=1,...,T
  w12[[1]] = w_0
  V12[[1]] = V_0 + diag(k)
  K = V12[[1]] %*% t(Beta[[1]]) %*% solve(Beta[[1]] %*% V12[[1]] %*% t(Beta[[1]]) + Sigma[[1]])
  w1[[1]] = w12[[1]] + K %*% (Y[,1] - Beta[[1]] %*% w12[[1]])
  V1[[1]] = V12[[1]] - K %*% Beta[[1]] %*% V12[[1]]
  
  for (t in 2:T){
    w12[[t]] = w1[[t-1]]
    V12[[t]] = V1[[t-1]] + diag(k)
    K = V12[[t]] %*% t(Beta[[t]]) %*% solve(Beta[[t]] %*% V12[[t]] %*% t(Beta[[t]]) + Sigma[[t]])
    w1[[t]] = w12[[t]] + K %*% (Y[,t] - Beta[[t]] %*% w12[[t]])
    V1[[t]] = V12[[t]] - K %*% Beta[[t]] %*% V12[[t]]
  }
 
# Backward Smoothing Step for t=1,...,T
  w[[T]] = w1[[T]]
  V[[T]] = V1[[T]]
  
  for (t in T:2){
    Z = V1[[t-1]] %*% solve(V12[[t]])
    w[[t-1]] = w1[[t-1]] + Z %*% (w[[t]] - w12[[t]])
    V[[t-1]] = V1[[t-1]] + Z %*% (V[[t]] - V12[[t]]) %*% t(Z)
  }
  
  
  result$w = w
  result$V = V
  result$w0 = w_0
  
  result
}

####################################################################

# Posterior probability of classifying beta into stationary slab

# Input: (all scalars)
#   beta: factor loading of interest (beta_{jk}^t)
#   lambda1, lambda0: slab & spike variance respectively
#   phi_beta: Autoregressive coefficient for beta
#   Theta: marginal importance weight

# Output:
#   posterior probability of classifying beta into stationary slab

#__________________________________________________________________


theta_beta = function(beta,lambda1,lambda0,phi_beta,Theta){
  
  num= dnorm(beta,0,sqrt(lambda1/(1-phi_beta^2)),log=T)
  
  den= dnorm(beta,0,sqrt(lambda0),log=T)
  
  1/(1+(1-Theta)/Theta*exp(den-num))
}


####################################################################

# Conditional inclusion probability of classifying beta into stationary slab, given previous value

# Input: (all scalars)
#   beta: factor loading of interest (beta_{jk}^t)
#   beta_previous: factor loading at previous time point (beta_{jk}^{t-1})
#   lambda1, lambda0: slab & spike variance respectively
#   phi_beta: Autoregressive coefficient for beta
#   theta: mixing weight (theta_{jk}^t)

# Output:
#   conditional inclusion probability into slab, given previous value

#__________________________________________________________________

pstar_beta = function(beta,beta_previous,lambda1,lambda0,phi_beta,theta){
  
  mu = phi_beta*beta_previous
  
  num= dnorm(beta,mu,sqrt(lambda1),log=T)
  
  den= dnorm(beta,0,sqrt(lambda0),log=T)
  
  1/(1+(1-theta)/theta*exp(den-num))
  
}

####################################################################

# Prior probability for Autoregressive coeffficient phi_beta

# Input: (all scalars)
#   phi_beta: value of AR coefficient for beta
#   a0, b0: parameters of a Beta(a0,b0) distribution

# Output:
#   prior probability for Autoregressive coeffficient phi_beta

#__________________________________________________________________

phi_beta_prior = function(phi_beta,a0,b0){
  
  phistar=(phi_beta+1)/2
  
  dbeta(phistar,a0,b0,log=TRUE)
  
}

####################################################################

# Log-likelihood for Gibbs sampling acceptance probability of phi_beta

# Input: (all scalars)
#   phi_beta: value of AR coefficient for beta
#   Beta: a list of T matrices; pxk matrix B_t = Beta[[t]]
#   Beta0: a pxk matrix B_0 (factor loading at t=0)
#   Gamma:  a list of T binary matrices; pxk matrix Gamma_t = Gamma[[t]]
#   Gamma0: a pxk matrix Gamma_0 (spike-slab indicator at t=0)
#   lambda1, lambda0: slab & spike variance respectively
#   Theta: marginal importance weight

# Output:
#   log-likelihood value

#__________________________________________________________________


loglik = function(phi_beta,Beta,Beta0,Gamma,Gamma0,lambda1,lambda0,Theta){
  
  p = nrow(Beta[[1]])
  T= length(Beta)
  value = 0
  
  for (j in 1:p){
    
    beta0 = Beta0[j,]
    gamma0 = Gamma0[j,]
    
    beta = Beta[[1]][j,]
    gamma = Gamma[[1]][j,]
    
    for (t in 2:T){
      beta = cbind(beta, Beta[[t]][j,])
      gamma = cbind(gamma, Gamma[[t]][j,])
    }
  
    thetas = theta_beta(beta0,lambda1,lambda0,phi_beta,Theta)
  
    thetas[thetas==1] = 1-10^{-10}
  
    thetas[thetas==0] = 10^{-10}
  
    value = value - sum(gamma0*beta0^2*(1-phi_beta^2)/(2*lambda1))
  
    value = value + 0.5*sum(gamma0)*log((1-phi_beta^2)/lambda1)
  
    value = value - sum(gamma[,1]*(beta[,1]-beta0*phi_beta)^2/(2*lambda1))
  
    value = value + sum(gamma[,1]*log(thetas)+(1-gamma[,1])*log(1-thetas))
  
    for (t in (2:T)){
    
      thetas = theta_beta(beta[,t-1],lambda1,lambda0,phi_beta,Theta)
    
      thetas[thetas==1] = 1-10^{-10}
    
      thetas[thetas==0] = 10^{-10}
    
      value = value - sum(gamma[,t]*(beta[,t]- beta[,t-1]*phi_beta)^2/(2*lambda1))
    
      value = value + sum(gamma[,t]*log(thetas) + (1-gamma[,t])*log(1-thetas))
    
    }
  }
  
  value
  
}

####################################################################

# Run MCMC for Dynamic Sparse Factor Analysis

# Input:
#   Y: pxT Data matrix
#   phi_w: autoregressive coefficient for latent factors
#   phi_beta_init: initial value of AR coefficient for beta
#   phiest: estimate phi_beta? (TRUE / FALSE) phi_beta = 1 if FALSE
#   lambda1, lambda0: slab & spike variance respectively
#   Theta: marginal importance weight
#   delta: discount factor for stochastic volatility
#   n0, d0: p-vectors, stochastic volatility parameters
#   RV: random walk on Beta? (TRUE / FALSE)
#   NITER: number of Monte Carlo iterations

# Output:
#   A list containing:
#     Sigma_sim: 
#         a list containing NITER-many MCMC draws of Sigma_{1:T}
#         Sigma is a list of T matrices; 
#         pxp matrix Sigma_t = Sigma_sim[[iter]][[t]]
#     Beta_sim: 
#         a list containing NITER-many MCMC draws of Beta_{1:T}
#         Beta is a list of T matrices; 
#         pxk matrix Beta_t = Beta_sim[[iter]][[t]]
#     Beta0_sim: 
#         a list containing NITER-many MCMC draws of Beta_0
#         Beta_0 is a pxk matrix
#     Gamma_sim:
#         a list containing NITER-many MCMC draws of Gamma_{1:T}
#         Gamma is a list of T matrices; 
#         pxk matrix Gamma_t = Gamma_sim[[iter]][[t]]
#     Gamma0_sim:
#         a list containing NITER-many MCMC draws of Gamma_0
#         Gamma_0 is a pxk matrix
#     w_sim:
#         a list containing NITER-many MCMC draws of latent factors w_{1:T}
#         k-vector w_t = w_sim[[iter]][[t]]
#     phi_beta_sim:
#         a list containing NITER-many MCMC draws of phi_beta
#         phi_beta is scalar
#     accept_sim: 
#         decision for Metropolis draws of  phi_beta (1 = accept, 0 = reject)
#     m, C: 
#         lists of p vectors & matrices respectively
#         posterior mean & variance for Beta
#         m[[j]], C[[j]] correspond to beta_{j.}^t
#     n, d : 
#         parameters for poseterior distributions of sigma_{jk}^t
#     RV: Random walk on factor loadings (phi_beta = 1)
#         TRUE / FALSE

#__________________________________________________________________


MCMC_DSS = function(Y, phi_w, phi_beta_init = 0.95, phiest = TRUE, lambda0,
                   lambda1, Theta, delta, n0, d0, RV = FALSE, NITER = 2000){
  
  p = nrow(Y)
  T = ncol(Y)
  k = floor(p*2/3)
  library(mvtnorm)
  
  # Initialization
  
  phi0 = 0
  phi_beta = phi_beta_init
  if (RV==TRUE) { phi_beta = 1 }
  
  Gamma0 = matrix(rbinom(1,1,Theta),p,k)
  Beta0 = matrix(rnorm(1,2,1), p, k)*Gamma0
  
  Gamma = list()
  Sigma = list()
  Beta = list()
  
  for (t in 1:T){ 
    Gamma[[t]] = matrix(rbinom(1,1,Theta),p,k)
    Sigma[[t]] = diag(p)
    Beta[[t]] = matrix(rnorm(1,2,1),p,k)*Gamma[[t]]
  }
  
  
  # Metropolis Hastings for phi_beta: 
  #   prior: Beta(a0,b0), proposal variance: var_proposal
  a0 = 20
  b0 = 1.5
  var_proposal = 0.001 
  
  # Sampled Factor Loadings
  Beta_sim = list()
  Beta0_sim = list()
  
  # Sampled spike-slab inclusion matrices
  Gamma_sim = list()
  Gamma0_sim = list()
  
  # Sampled dynamic covariance matrices
  Sigma_sim = list()
  
  # Sampled latent factors
  w_sim = list()
  
  # Sampled AR coefficients: phi_beta
  phi_beta_sim = list()
  
  # (1 = accept, 0 = reject) decision for Metropolis draws of  phi_beta
  accept_sim = list()
  
  Y_est = list()
  for (iter in 1:NITER){Y_est[[iter]] = matrix(0,p,T)}
  
  NAS = FALSE
  
  # Start MCMC iterations:
  iter = 0
  X = list()
  
  while (iter < NITER & NAS == FALSE){
    
    iter = iter + 1
    cat("Iteration: ", iter, "\n")
    
    # Kalman filter & smoother for latent factors
    KF_result = Kalman_FS(phi_w, Beta, Sigma, Y)
    for (t in 1:T){
      CH = chol(KF_result$V[[t]], pivot=TRUE)
      #X[[t]] = KF_result$w[[t]] + t(CH) %*% t(rmvnorm(1, rep(0,k), diag(k)))
      X[[t]] = solve(t(CH)) %*% KF_result$w[[t]] + t(rmvnorm(1, rep(0,k), diag(k)))
    }
    
    
    
    ########## Sampling Rotation Matrix ##########
    
    #X0 = KF_result$w0
    #V_rotate = matrix(0,k,k)
    #(X[[1]] - X0) %*% t((X[[1]] - X0))
    #for (t in 2:T){
     # V_rotate = V_rotate + (X[[t]] - X[[t-1]]) %*% t((X[[t]] - X[[t-1]]))
    #}
    #V_rotate = V_rotate/T
    #V_L = chol(V_rotate, pivot=TRUE)
    #A = V_L %*% rWishart(1, k, diag(k))[,,1] %*% t(V_L)
    #A_L = chol(A, pivot=TRUE)
    
    ########## Sampling Factor Loadings: Beta_{1:T} ##########
    
    # Forward Filtering
    
    a = list()
    R = list()
    C = list()
    m = list()
    e = list()
    q = list()
    f = list()
    phi = list()
    for (t in 1:T){
      a[[t]] = list()
      R[[t]] = list()
      m[[t]] = list()
      C[[t]] = list()
      e[[t]] = list()
      q[[t]] = list()
      f[[t]] = list()
    }
    m0 = list()
    C0 = list()
    
    for (j in 1:p){
      m0[[j]] = phi0 * Gamma0[j,]
      if (RV == FALSE){
        C0[[j]] = diag(Gamma0[j,]*lambda1/(1-phi_beta^2)+(1-Gamma0[j,])*lambda0)
      } else{
        C0[[j]] = diag(Gamma0[j,]*lambda1+(1-Gamma0[j,])*lambda0)}
    }
    
    for (t in (1:T)){
      for (j in 1:p){
      
      H = Gamma[[t]][j,] * phi0
      W = diag(Gamma[[t]][j,] * lambda1 + (1-Gamma[[t]][j,])*lambda0)
      G = diag(Gamma[[t]][j,] * phi_beta)
      
      if (t > 1){
        a[[t]][[j]] = H + G %*% (m[[t-1]][[j]]-H)
        R[[t]][[j]] = G %*% C[[t-1]][[j]] %*% t(G) + W
      } else{
        a[[t]][[j]] = H + G %*% (m0[[j]]-H)
        R[[t]][[j]] = G %*% C0[[j]] %*% t(G) + W
      }
      
      f[[t]][[j]] = t(a[[t]][[j]]) %*% X[[t]]
      e[[t]][[j]] = as.numeric(Y[j,t]-f[[t]][[j]])
      
      q[[t]][[j]] = as.numeric(t(X[[t]]) %*% R[[t]][[j]] %*% X[[t]]) + Sigma[[t]][j,j]	
      A = R[[t]][[j]] %*% X[[t]] / q[[t]][[j]]
      m[[t]][[j]] = a[[t]][[j]] + e[[t]][[j]] * A
      C[[t]][[j]] = R[[t]][[j]] - A  %*% t(A) * q[[t]][[j]]
      }
    }
    
    # Backwards Smoothing
    
    for (t in (T:1)){
      B = list()
      for (j in 1:p){
        if (t<T){
      
          B = C[[t]][[j]] %*% diag(Gamma[[t+1]][j,]*phi_beta) %*% solve(R[[t+1]][[j]])
          mean = m[[t]][[j]] + B %*% ( Beta[[t+1]][j,] - a[[t+1]][[j]])
          variance = C[[t]][[j]] - B %*% R[[t+1]][[j]] %*% t(B)
        } else{
        
        mean = m[[t]][[j]]
        variance = C[[t]][[j]]
      }
        CH = chol(variance, pivot=TRUE)	
        Beta[[t]][j,] = mean + t(CH) %*% rnorm(k)
      }
    }
    
    
    ########## Sample Beta0 ##########
    
    for (j in 1:p){
      
      B = C0[[j]] %*% diag(Gamma[[1]][j,] * phi_beta) %*% solve(R[[1]][[j]])
      mean = m0[[j]] + B %*% ( Beta[[1]][j,] - a[[1]][[j]])
      variance = C0[[j]] - B %*% R[[1]][[j]] %*% t(B)
      CH = chol(variance, pivot=TRUE)	
      Beta0[j,] = mean + t(CH) %*% rnorm(k) 
      
    }
    
    
    #for (i in 1:(k_max-1)){for (j in (i+1):k_max){Beta0[i,j]=0}}
    #for (t in 1:T){
     # for (i in 1:(k_max-1)){for (j in (i+1):k_max){Beta[[t]][i,j]=0}}
    #}
    
    Beta_sim[[iter]] = Beta
    Beta0_sim[[iter]] = Beta0
    
    w_sim[[iter]] = X
    
    ########## Sample Gamma_{1:T} ##########
    
    for (t  in (1:T)){
      for (j in 1:p){
        
        if (t>1){
          if (RV==FALSE){
            thetas = theta_beta(Beta[[t-1]][j,],
                               lambda1, lambda0, phi_beta, Theta)
            ps = pstar_beta(Beta[[t]][j,], Beta[[t-1]][j,],
                           lambda1, lambda0, phi_beta, thetas)
            }
          else{ 
            ps = pstar_beta(Beta[[t]][j,], Beta[[t-1]][j,],
                            lambda1, lambda0, phi_beta, Theta)
          }
        }
        else{
          if (RV==FALSE){
            thetas = theta_beta(Beta0[j,],
                               lambda1, lambda0, phi_beta, Theta)
            ps = pstar_beta(Beta[[t]][j,], Beta0[j,],
                           lambda1, lambda0, phi_beta, thetas)
            }
          else{
            ps = pstar_beta(Beta[[t]][j,], Beta0[j,],
                           lambda1, lambda0, phi_beta, Theta)
            }
        }
        Gamma[[t]][j,] = rbinom(k,1,ps)
      }
    }
    
    ########## Sample Gamma0 ##########
    
    for (j in 1:p){
      
      if (RV==FALSE){
        thetas = theta_beta(Beta0[j,], lambda1, lambda0, phi_beta, Theta)
      }
      else{
        thetas = pstar_beta(Beta0[j,], numeric(p), lambda1, lambda0, phi_beta, Theta)
      }
      
      Gamma0[j,] = rbinom(k,1,thetas)
    }
      
    Gamma_sim[[iter]] = Gamma
    Gamma0_sim[[iter]] = Gamma0
      
    ########## Sample variances Sigma_{1:T} ##########
      
    n = list()
    d = list()
    phi = list()
    
    for (t in 1:T){
      n[[t]] = numeric(p)
      d[[t]] = numeric(p)
      phi[[t]] = numeric(p)
      }

    # Forward filtering
    
    for (j in 1:p){
      for (t in (1:T)){
        
        beta_th = Beta[[t]][j,]
        
        if (t>1){
          r =  as.numeric(Y[j,t]- t(X[[t]]) %*% beta_th)
          n[[t]][j] = delta * n[[t-1]][j] + 1
          d[[t]][j] = delta * d[[t-1]][j] + r^2
          }
        else{
          r = as.numeric(Y[j,t]- t(X[[t]]) %*% beta_th)
          n[[t]][j] = delta * n0[j] + 1
          d[[t]][j] = delta * d0[j] + r^2
        }
      }
      
    # Backwards Smoothing
      
      for (t in (T:1)){
        
        if (t==T){
          shape = n[[T]][j]/2
          rate = d[[T]][j]/2
          phi[[T]][j] = rgamma(1,shape,rate)
        }
        if (t<T){
          shape = (1-delta) * n[[t]][j]/2
          rate = d[[t]][j]/2
          nu = rgamma(1,shape=shape,rate=rate)
          phi[[t]][j] = nu + delta * phi[[t+1]][j]
        }
      }
    }
    
    for (t in 1:T){
      Sigma[[t]] = diag(1/phi[[t]])
    }
    
    Sigma_sim[[iter]] = Sigma
    
    ########## Sample phi_beta ##########
    
    if ( (RV==FALSE) & (phiest==TRUE)){
      
      phi_beta_new = 0.95 + runif(1)*0.05
      
      if (phi_beta_new > 0.999){phi_beta_new = 0.999}
      
      log_ratio = loglik(phi_beta_new, Beta, Beta0, Gamma,
                          Gamma0, lambda1, lambda0, Theta) +
                  phi_beta_prior(phi_beta_new, a0, b0) -
                  loglik(phi_beta, Beta, Beta0, Gamma,
                          Gamma0, lambda1, lambda0, Theta) -
                  phi_beta_prior(phi_beta, a0, b0)
      
      alpha = min(1, exp(log_ratio))
      accept = rbinom(1, 1, alpha)
      phi_beta = phi_beta_new * accept + (1-accept) * phi_beta
      accept_sim[[iter]] = accept
    }
    phi_beta_sim[[iter]] = phi_beta
    
    ########## Estimate Y ##########
    
    for (t in 1:T){
      Y_est[[iter]][,t] = Beta_sim[[iter]][[t]] %*% w_sim[[iter]][[t]] 
    }
    
    nas = 0
    for (t in 1:T){
      nas = nas + sum(apply(Beta[[t]], 2, is.na))
    }
    NAS = nas > 0
  }
  list(Sigma_sim = Sigma_sim,
       Beta_sim = Beta_sim,
       Beta0_sim = Beta0_sim,
       Gamma_sim = Gamma_sim,
       Gamma0_sim = Gamma0_sim,
       w_sim = w_sim,
       phi_beta_sim = phi_beta_sim,
       accept_sim = accept_sim,
       m = m[[T]],
       C = C[[T]],
       n = n,
       d = d,
       RV = RV,
       Y_est = Y_est)
}


####################################################################

# Compute mean from MCMC draws

# Input:
#   B: a list of MCMC draws, 
#      each draw is a list of T elements (eg. Beta_{1:T}, Sigma_{1:T}, phi_beta etc)
#   burin: number of draws to be removed as burn-in period

# Output:
#   Posterior mean estimated from MCMC

#__________________________________________________________________

Parse = function(B, burnin=500){
  
  B_type = class(B[[1]][[1]])
  T = length(B[[1]])
  
  if (B_type == "matrix"){
    p1 <- nrow(B[[1]][[1]])
    p2 <- ncol(B[[1]][[1]])
    
    Baggreg = list()
    
    for (t in 1:T){
      Baggreg[[t]] = matrix(0, p1, p2)
      for (i in burnin:length(B)){
        Baggreg[[t]] = Baggreg[[t]] + B[[i]][[t]]
      }
      Baggreg[[t]] = Baggreg[[t]]/(length(B)-burnin)
    }
  } else{
    
    B_type = class(B[[1]])
    
    if (B_type =="matrix"){
      p1 <- nrow(B[[1]])
      p2 <- ncol(B[[1]])
      Baggreg = matrix(0, p1, p2)
      for (i in burnin:length(B)){
        Baggreg = Baggreg + B[[i]]
      }
      Baggreg = Baggreg/(length(B)-burnin)
    } else{
      B_type = class(B[[1]])
      
      if (B_type =="matrix"){
        p1 <- nrow(B[[1]])
        p2 <- ncol(B[[1]])
        Baggreg = matrix(0, p1, p2)
        for (i in burnin:length(B)){
          Baggreg = Baggreg + B[[i]]
        }
        Baggreg = Baggreg/(length(B)-burnin)
      }
      Baggreg = 0
      for (i in (burnin:length(B))){
        Baggreg = Baggreg + B[[i]]
      }
      Baggreg = Baggreg/(length(B)-burnin)
    }
  }
  
  Baggreg
}

####################################################################

# Find credible region from MCMC draws

# Input:
#   B: a list of MCMC draws, 
#      each draw is a list of T elements (eg. Beta_{1:T}, Sigma_{1:T}, phi_beta etc)
#   burin: number of draws to be removed as burn-in period

# Output: a list containing
#     upper 25% quantile
#     lower 25% quantile
#     median

#__________________________________________________________________

Credible = function(B,burnin=500){
  
  quantile_u = B[[1]]
  quantile_l = B[[1]]
  median = B[[1]]
  
  B_type = class(B[[1]][[1]])
  
  if (B_type == "matrix"){
    p1 <- nrow(B[[1]][[1]])
    p2 <- ncol(B[[1]][[1]])
    for (t in 1:T){
      for (i in 1:p1){
        for (j in 1:p2){
          mcmc_sample = vector()
          for (iter in burnin:length(B)){
            mcmc_sample = append(mcmc_sample, B[[iter]][[t]][i,j])
          }
          quantile_u[[t]][i,j] = as.numeric(quantile(mcmc_sample,0.975))
          quantile_l[[t]][i,j] = as.numeric(quantile(mcmc_sample,0.025))
          median[[t]][i,j] = as.numeric(quantile(mcmc_sample,0.5))
        }
      }
    }
  } else{
    
    B_type = class(B[[1]])
    
    if (B_type =="matrix"){
      p1 <- nrow(B[[1]])
      p2 <- ncol(B[[1]])
      for (i in 1:p1){
        for (j in 1:p2){
          mcmc_sample = vector()
          for (iter in burnin:length(B)){
            mcmc_sample = append(mcmc_sample, B[[iter]][i,j])
          }
          quantile_u[i,j] = as.numeric(quantile(mcmc_sample,0.975))
          quantile_l[i,j] = as.numeric(quantile(mcmc_sample,0.025))
          median[i,j] = as.numeric(quantile(mcmc_sample,0.5))
        }
      }
    } else{
      mcmc_sample = vector()
      for (iter in burnin:length(B)){
        mcmc_sample = append(mcmc_sample, B[[iter]])
      }
      quantile_u = as.numeric(quantile(mcmc_sample,0.975))
      quantile_l = as.numeric(quantile(mcmc_sample,0.025))
      median = as.numeric(quantile(mcmc_sample,0.5))
    }
  }
  
  list(upper = quantile_u,
       lower = quantile_l,
       median = median)
  
}

####################################################################

# Threshold Beta to 0, if less than 'thres' 

# Input:
#   B: a list of T matrices Beta_{1:T}
#   thres: spike-threshold, usually set to lambda0 (spike sd)

# Output: a list of T matrices Beta_{1:T}, sparse

#__________________________________________________________________

sparsify_beta = function(B, thres = 0.1){
  T = length(B)
  p = nrow(B[[1]])
  k = ncol(B[[1]])
  for (t in 1:T){
    B[[t]] = B[[t]] * (B[[t]] > thres)
  }
  B
}

####################################################################

# One-step-ahead Forecast

# Input:
#   Y: pxT data metrix
#   result: MCMC_DFA output, a list with simulated Beta, Sigma, w, phi_beta
#   lambda0,lambda1: spike & slab variances
#   Theta: marginal importance weight
#   delta: stochastic volatility parameter
#   burnin: number of MCMC draws to be removed as burn-in period
#   RV: random walk on Beta? (TRUE/FALSE)

# Output: One-step-ahead forecast, a p-vector Y_{T+1}

#__________________________________________________________________

forecast_MCMC = function(Y, results, lambda1, lambda0, phi_w = 0.98, delta, Theta, burnin, RV){
  
  if(RV==FALSE){ phihat = mean(unlist(results$phi_beta_sim)[-(1:burnin)])} else {phihat = 1}
  
  mean = results$m
  variance = results$C
  
  Bmeans = Parse(results$Beta_sim, burnin)
  Sigma_means = Parse(results$Sigma_sim, burnin)
  
  T = length(results$d)
  p = nrow(Bmeans[[1]])
  k = ncol(Bmeans[[1]])
  
  n0 = rep(1,p)
  d0 = n0
  
  X = Parse(results$w_sim, burnin)
  
  n = list()
  d = list()
  for (t in 1:T){
    n[[t]] = numeric(p)
    d[[t]] = numeric(p)
  }
  
  # Forward filtering
  
  xnew = phi_w * X[[T]]
  #Vnew = diag(k)
  #K = Vnew %*% t(Bmeans[[T]]) %*% solve(Bmeans[[T]] %*% Vnew %*% t(Bmeans[[T]]) + Sigma_means[[T]])
  #xnew = xnew + K %*% (Y[,T] - Bmeans[[T]] %*% xnew)
  
  for (t in (1:T)){
    for (j in 1:p){
      beta_th = Bmeans[[t]][j,]
      if (t>1){
        r =  as.numeric(Y[j,t]- t(X[[t]]) %*% beta_th)
        n[[t]][j] = delta * n[[t-1]][j] + 1
        d[[t]][j] = delta * d[[t-1]][j] + r^2
      }else{
      r =  as.numeric(Y[j,t]- t(X[[t]]) %*% beta_th)
      n[[t]][j] = delta * n0[j] + 1
      d[[t]][j] = delta * d0[j] + r^2
      }
    }
  }
  
  gamma_new = list()
  a_new = list()
  q = numeric(p)
  f_new = numeric(p)
  f2_new = numeric(p)
  
  for (j in 1:p){
    if(RV==FALSE){
      gamma_new[[j]] = theta_beta(mean[[j]], lambda1, lambda0, phihat, Theta)
    } else{
      gamma_new = Theta
    }
    a_new[[j]] = gamma_new[[j]] * phihat * mean[[j]]
    q[j] = as.numeric(t(xnew) %*% variance[[j]] %*% xnew + d[[T]][j]/n[[T]][j])
    f_new[j] = as.numeric(t(xnew) %*% a_new[[j]])
    f2_new[j] = as.numeric(t(xnew) %*% mean[[j]])
  }
  
  if (RV==FALSE){pred_mean = f_new} else {pred_mean = f2_new}
  
  list(type1 = f_new, type2 = f2_new, pred_mean = pred_mean, pred_variance = q)
  
}
