# Dynamic Sparse Factor Analysis: EM Implementation
# Call:
#   results = EM_DSS(Y, phi_w, phi_beta = 0.95,
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

Kalman_FS = function(phi_w, Beta, Sigma, Y){
  
  T = length(Beta)
  p = nrow(Beta[[1]])
  k = ncol(Beta[[1]])
  
  # Initialize w_0|0 and V_0|0
  w_0 = rep(0,k)
  V_0 = 1/(1-phi_w^2)*diag(k)
  
  for (iter in 1:20){
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
    K = V12[[1]] %*% t(Beta[[1]]) %*% solve(Beta[[1]] %*% V12[[1]] %*% t(Beta[[1]]) + Sigma[[1]])
    w1[[1]] = w12[[1]] + K %*% (Y[,1] - Beta[[1]] %*% w12[[1]])
    V1[[1]] = V12[[1]] - K %*% Beta[[1]] %*% V12[[1]]
    V1[[1]] = (V1[[1]] + t(V1[[1]]))/2
    
    for (t in 2:T){
      w12[[t]] = w1[[t-1]]
      V12[[t]] = V1[[t-1]] + diag(k)
      K = V12[[t]] %*% t(Beta[[t]]) %*% solve(Beta[[t]] %*% V12[[t]] %*% t(Beta[[t]]) + Sigma[[t]])
      w1[[t]] = w12[[t]] + K %*% (Y[,t] - Beta[[t]] %*% w12[[t]])
      V1[[t]] = V12[[t]] - K %*% Beta[[t]] %*% V12[[t]]
      V1[[t]] = (V1[[t]] + t(V1[[t]]))/2
    }
    
    # Backward Smoothing Step for t=1,...,T
    w[[T]] = w1[[T]]
    V[[T]] = V1[[T]]
    Cov[[T]] = (diag(k) - K %*% Beta[[T]]) %*% V1[[T]]
    
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
  }
  
  result$w = w
  result$V = V
  result$w0 = w_0
  result$V0 = V_0
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

####################################################################

# Log-likelihood

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

# Run EM for Dynamic Sparse Factor Analysis

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

#__________________________________________________________________

EM_DSS = function(Y, phi_w, phi_beta = 0.95,
                  lambda0, lambda1, Theta, delta, n0, d0, 
                  NITER = 10000, epsilon = 0.1, PX = FALSE){
  
  
  p = nrow(Y)
  T = ncol(Y)
  k = floor(p*2/3)
  library(mvtnorm)
  
  # Initialization
  
  phi0 = 0
  
  Gamma0 = matrix(rbinom(p*k,1,Theta),p,k)
  Beta0 = matrix(rnorm(p*k,2,1), p, k)*Gamma0
  
  Gamma = list()
  Sigma = list()
  Beta = list()
  
  for (t in 1:T){ 
    Gamma[[t]] = matrix(rbinom(p*k,1,Theta),p,k)
    Sigma[[t]] = diag(p)
    Beta[[t]] = matrix(rnorm(p*k,2,1),p,k)*Gamma[[t]]
  }
  
  #Beta = Beta_true
  #Beta0 = Beta0_true
  
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
  Beta0_new = Beta0
  
  eps = epsilon+1
  
  while((eps > epsilon) &  (iter < NITER)){
    
    iter = iter + 1
    cat("Iteration: ", iter, "\n")
    
    ####*****###### The E-Step ####*****######
    
    ##### Update Latent Factors: Kalman Filter #####
    
    KF_result = Kalman_FS(phi_w=0.98, Beta, Sigma, Y)
    w = KF_result$w
    
    ##### Update Indicators: Gamma_{1:T} #####
    
    for (t  in (1:T)){
      for (j in 1:p){
        
        if (t>1){
          thetas = theta_beta(Beta[[t-1]][j,],
                              lambda1, lambda0, phi_beta, Theta)
          ps = pstar_beta(Beta[[t]][j,], Beta[[t-1]][j,],
                          lambda1, lambda0, phi_beta, thetas)
        }else{
          thetas = theta_beta(Beta0[j,],
                              lambda1, lambda0, phi_beta, Theta)
          ps = pstar_beta(Beta[[t]][j,], Beta0[j,],
                          lambda1, lambda0, phi_beta, thetas)
        }
        Gamma[[t]][j,] = ps
      }
    }
    
    ########## Update Gamma0 ##########
    
    for (j in 1:p){
      thetas = theta_beta(Beta0[j,], lambda1, lambda0, phi_beta, Theta)
      Gamma0[j,] = thetas
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
        
        beta_th = Beta[[t]][j,]
        beta_th[Gamma[[t]][j,] < 0.5]= 0
        
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
    
    #w_tilde = list()
    
    #for (t in 1:T){
    #  svd_latent = svd(KF_result$V[[t]])
    #  U_tilde = diag(sqrt(svd_latent$d)) %*% svd_latent$u
    #  w_tilde[[t]] = matrix(0,k+1,k)
    #  w_tilde[[t]][1,] = KF_result$w[[t]]
    #  w_tilde[[t]][2:(k+1),] = U_tilde
    #}
    
    for (j in 1:p){
      S = list()
      V = list()
    
      for(t in (1:T)){
        V[[t]] = Gamma[[t]][j,]/lambda1+(1-Gamma[[t]][j,])/lambda0
        D = V[[t]]		
        if(t<T){
          D = D + phi_beta^2*(Gamma[[t+1]][j,])/lambda1 
        }
        xrow = w[[t]]/D
        S[[t]] = diag(1/D)-phi[[t]][j]*xrow%*%t(xrow)/(1+phi[[t]][j]*sum(w[[t]]^2/D))
      }
  
      nrep = 100
      for (loop in (1:nrep)){
        permute = sample(1:T)
        
        for (t in permute){
          if (t>1){ previous = Beta[[t-1]][j,] } else { previous=Beta0[j,] } 
          mean = phi[[t]][j]*w[[t]]*Y[j,t] + phi_beta*previous*Gamma[[t]][j,]/lambda1
        
          if(t<T){
            mean = mean + phi_beta*Beta[[t+1]][j,]*Gamma[[t+1]][j,]/lambda1
          }
        
          Beta[[t]][j,] = S[[t]] %*% mean
        
        }
        # Update Beta0
      
        num = Gamma[[1]][j,] * Beta[[1]][j,] * phi_beta / lambda1
        denom = Gamma0[j,]/(lambda1/(1-phi_beta^2)) +
          (1-Gamma0[j,])/lambda0 + phi_beta^2 * Gamma[[1]][j,]/lambda1
        Beta0[j,] = num/denom
      }
    }
    
    ########## Update Rotation Matrices A_{1:T} ##########
    if (PX == TRUE){
      A = list()
      M1 = KF_result$w0 %*% t(KF_result$w0) + KF_result$V0
      M12 = KF_result$w0 %*% t(KF_result$w[[1]]) + KF_result$Cov[[1]]
      M21 = t(M12)
      M2 = KF_result$w[[1]] %*% t(KF_result$w[[1]]) + KF_result$V[[1]]
      A[[1]] = M1 - M12 - M21 + M2
      
      for (t in 2:T){
        M1 = KF_result$w[[t-1]] %*% t(KF_result$w[[t-1]]) + KF_result$V[[t-1]]
        M12 = KF_result$w[[t-1]] %*% t(KF_result$w[[t]]) + KF_result$Cov[[t]]
        M21 = t(M12)
        M2 = KF_result$w[[t]] %*% t(KF_result$w[[t]]) + KF_result$V[[t]]
        A[[t]] = M1 - M12 - M21 + M2
      }
      
      
      ########## Rotation to Sparsity ##########
      
      A_L = list()
      for (t in 1:T){
        A_L[[t]] = chol(A[[t]], pivot=TRUE)
        Beta[[t]] = Beta[[t]] %*% A_L[[t]]
      }
    }
    
    
    ########## Evaluate stopping criteria ##########
    
    eps = max(abs(Beta0_new-Beta0))
    for (t in 1:T){
      eps = max(eps,max(abs(Beta_new[[t]]-Beta[[t]])))
    }
    
    
    Beta_new = Beta	
    Beta0_new = Beta0
    
    print(eps)
  }
  
  list(Sigma = Sigma,
       Beta = Beta,
       Beta0 = Beta0,
       Gamma = Gamma,
       Gamma0 = Gamma0,
       w = w)
}


