# Dynamic Sparse Factor Analysis: EM Implementation (Static Loadings)
# Call:
#   results = EM_static(Y, phi_w, phi_beta = 0.95,
#     lambda0, lambda1, Theta, delta, n0, d0, 
#     NITER = 10000, epsilon = 0.1)

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
#       result$Cov: a list of T matrices; kxk matrix V_{t,t-1 | T} = result$Cov[[t]]


# Other Notations:
#   w1: a list of T vectors: k-vector w_{t|t} = w1[[t]]
#   w12: a list of T vectors: k-vector w_{t|t-1} = w12[[t]]
#   V1: a list of T matrices: kxk matrix V_{t|t} = V1[[t]]
#   V12: a list of T matrices: kxk matrix V_{t|t-1} = V12[[t]]

#__________________________________________________________________

Kalman_FS_static = function(phi_w, Beta, Sigma, Y){
  
  iter = 10
  
  T = length(Sigma)
  p = nrow(Beta)
  k = ncol(Beta)
  
  # Initialize w_0|0 and V_0|0
  w_0 = rep(0,k)
  V_0 = 1/(1-phi_w^2)*diag(k)
  
  for (iter in 1:iter){
    w1 = list()
    w12 = list()
    V1 = list()
    V12 = list()
    
    result = list()
    w = list()
    V = list()
    Cov = list()
    
    # Forward Prediction & Correction Step for t=1,...,T
    w12[[1]] = w_0 
    V12[[1]] = V_0 + diag(k)
    K = V12[[1]] %*% t(Beta) %*% solve(Beta %*% V12[[1]] %*% t(Beta) + Sigma[[1]])
    w1[[1]] = w12[[1]] + K %*% (Y[,1] - Beta %*% w12[[1]])
    V1[[1]] = V12[[1]] - K %*% Beta %*% V12[[1]]
    V1[[1]] = (V1[[1]] + t(V1[[1]]))/2
    
    for (t in 2:T){
      w12[[t]] = w1[[t-1]]
      V12[[t]] = V1[[t-1]] + diag(k)
      K = V12[[t]] %*% t(Beta) %*% solve(Beta %*% V12[[t]] %*% t(Beta) + Sigma[[t]])
      w1[[t]] = w12[[t]] + K %*% (Y[,t] - Beta %*% w12[[t]])
      V1[[t]] = V12[[t]] - K %*% Beta %*% V12[[t]]
      V1[[t]] = (V1[[t]] + t(V1[[t]]))/2
    }
    
    # Backward Smoothing Step for t=1,...,T
    w[[T]] = w1[[T]]
    V[[T]] = V1[[T]]
    Cov[[T]] = (diag(k) - K %*% Beta) %*% V1[[T]]
    
    for (t in T:3){
      Z = V1[[t-1]] %*% solve(V12[[t]])
      w[[t-1]] = w1[[t-1]] + Z %*% (w[[t]] - w12[[t]])
      V[[t-1]] = V1[[t-1]] + Z %*% (V[[t]] - V12[[t]]) %*% t(Z)
      Z1 = V1[[t-2]] %*% solve(V12[[t-1]]) * phi_w
      Cov[[t-1]] = V1[[t-1]] %*% t(Z1) + Z %*% (Cov[[t]] - V1[[t-1]]) %*% t(Z1)
    }
    Z = V1[[1]] %*% solve(V12[[2]]) * phi_w
    w[[1]] = w1[[1]] + Z %*% (w[[2]] - w12[[2]])
    V[[1]] = V1[[1]] + Z %*% (V[[2]] - V12[[2]]) %*% t(Z)
    Z1 = V_0 %*% solve(V12[[1]]) * phi_w
    Cov[[1]] = V1[[1]] %*% t(Z1) + Z %*% (Cov[[2]] - V1[[1]]) %*% t(Z1)
    
    Z = V_0 %*% solve(V12[[1]]) * phi_w
    w_0 = w_0 + Z %*% (w[[1]] - w12[[1]])
    V_0 = V_0 + Z %*% (V[[1]] - V12[[1]]) %*% t(Z)
    V_0 = (V_0 + t(V_0))/2
  }
  
  result$w = w
  result$V = V
  result$w_0 = w_0
  result$V_0 = V_0
  result$Cov = Cov
  
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

# Run EM for Static Sparse Factor Analysis

# Input:
#   Y: pxT Data matrix
#   phi_w: autoregressive coefficient for latent factors
#   phi_beta: value of AR coefficient for beta
#   lambda1, lambda0: slab & spike variance respectively
#   Theta: marginal importance weight
#   delta: discount factor for stochastic volatility
#   n0, d0: p-vectors, stochastic volatility parameters
#   NITER: number of maximum iterations
#   epsilon: threshold for stopping EM
#   PX = TRUE/ FALSE: Parameter expanded model
#   const_rot = TRUE/ FALSE: rotation matrix constant over time or not

# Output:
#   A list containing:
#     Sigma: 
#         a list containing pxp matrix Sigma_t = Sigma[[t]]
#     Beta: 
#         a list containing pxk matrix Beta_t = Beta[[t]]
#     Beta0: 
#         Beta0 is a pxk matrix
#     Gamma:
#         a list containing pxk matrix Gamma_t = Gamma[[t]]
#     Gamma0:
#         Gamma0 is a pxk matrix
#     w:
#         a list containing latent factors k-vector w_t = w[[t]]
#     

#__________________________________________________________________

EM_static = function(Y, phi_w, phi_beta = 0.95,
                  lambda0, lambda1, Theta, delta, n0, d0, 
                  NITER = 10000, epsilon = 0.1){
  
  
  p = nrow(Y)
  T = ncol(Y)
  k = floor(p*2/3)
  library(mvtnorm)
  
  # Initialization
  
  phi0 = 0
  
  Gamma = matrix(rbinom(p*k,1,Theta),p,k)
  Beta = matrix(rnorm(p*k,2,1), p, k)*Gamma0
  
  Sigma = list()
  
  for (t in 1:T){ 
    Sigma[[t]] = diag(p)
  }
  
  w_0 = rep(0,k)
  V_0 = 1/(1-phi_w^2)*diag(k)
  
  # Estimated latent factors
  w = list()
  
  # Augmented data for M-step
  Y_tilde = list()
  for (t in 1:T){
    Y_tilde[[t]] = matrix(0, k+1, p)
    Y_tilde[[t]][1,] = Y[,t]
  }
  
  iter = 0
  
  Beta_new = Beta
  
  eps = epsilon+1
  
  while((eps > epsilon) &  (iter < NITER)){
    
    iter = iter + 1
    cat("Iteration: ", iter, "\n")
    
    ####*****###### The E-Step ####*****######
    
    ##### Update Latent Factors: Kalman Filter #####
    
    KF_result = Kalman_FS_static(phi_w=0.98, Beta, Sigma, Y)
    w = KF_result$w
    w_0 = KF_result$w_0
    V_0 = KF_result$V_0
    
    
    ##### Update Indicators: Gamma #####
    
    for (j in 1:p){
      thetas = theta_beta(Beta[j,], lambda1, lambda0, phi_beta, Theta)
      Gamma[j,] = thetas
    }
    ########## Update Sigma_{1:T} ##########
    
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
        
        beta_th = Beta[j,]
        beta_th[Gamma[j,] < 0.5] = 0
        
        if (t>1){
          r =  as.numeric(Y[j,t]- t(KF_result$w[[t]]) %*% beta_th)
          n[[t]][j] = delta * n[[t-1]][j] + 1
          d[[t]][j] = delta * d[[t-1]][j] + r^2
        }
        else{
          r = as.numeric(Y[j,t]- t(KF_result$w[[t]]) %*% beta_th)
          n[[t]][j] = delta * n0[j] + 1
          d[[t]][j] = delta * d0[j] + r^2
        }
      }
      
      # Backwards Smoothing
      
      for (t in (T:1)){
        
        if (t==T){
          shape = n[[T]][j]/2
          rate = d[[T]][j]/2
          phi[[T]][j] = shape/rate
        }
        if (t<T){
          shape = (1-delta) * n[[t]][j]/2
          rate = d[[t]][j]/2
          nu = shape/rate
          phi[[t]][j] = nu + delta * phi[[t+1]][j]
        }
      }
    }
    
    for (t in 1:T){
      Sigma[[t]] = diag(1/phi[[t]])
    }
    
    ####*****###### The M-Step ####*****######
    
    ########## Sampling Factor Loadings: Beta_{1:T} ##########
    
    w_tilde = list()
    
    for (t in 1:T){
      svd_latent = svd(KF_result$V[[t]])
      U_tilde = diag(sqrt(svd_latent$d)) %*% svd_latent$u
      w_tilde[[t]] = matrix(0,k+1,k)
      w_tilde[[t]][1,] = KF_result$w[[t]]
      w_tilde[[t]][2:(k+1),] = U_tilde
    }
    
    for (j in 1:p){
      V = list()
      mean = list()
      D = Gamma[j,]/lambda1+(1-Gamma[j,])/lambda0
      for(t in (1:T)){
        xrow = t(w_tilde[[t]])
        V[[t]] = diag(D) + xrow %*% t(xrow) * phi[[t]][j]
      }
      
      nrep = 1
      
        for (t in 1:T){
          mean[[t]] = phi[[t]][j] * w[[t]] * Y[j,t]
        }
        
        mn = mean[[t]]*0
        v = V[[t]]*0
        
        for (t in 1:T){
          mn = mn + mean[[t]]
          v = v + V[[t]]
        }
        
        Beta[j,] = solve(v) %*% mn
      }
    
    ########## Evaluate stopping criteria ##########

    eps = max(abs(Beta_new-Beta))
    
    Beta_new = Beta
    
    print(eps)
  }
  
  list(Sigma = Sigma,
       Beta = Beta,
       Gamma = Gamma,
       w = w)
}

####################################################################

# One-step-ahead Forecast

# Input:
#   result: EM_static output, a list with simulated Beta, Sigma, w, phi_beta
#   phi_w: AR coefficient of latent factors

# Output: One-step-ahead forecast, a p-vector Y_{T+1}

#__________________________________________________________________

forecast_SEM = function(result,phi_w){
  
  Beta = result$Beta
  w = result$w
  
  T = ncol(Beta)
  w_current = w[[T]]
  
  w_new = phi_w*w_current
  Y_new = Beta %*% w_new
  
  Y_new
  
}
