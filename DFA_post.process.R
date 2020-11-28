####################################################################

# Post processing for factor identification (Kaufman 2019)

# Input:
#   B: a list of MCMC draws of factor loadings, 
#   burin: number of draws to be removed as burn-in period

# Output:
#   Posterior mean estimated from MCMC after post processing

#__________________________________________________________________


post_process = function(B, w, burnin){
  
  # Set post processing parameters
  b_corr = 0.7
  b_draws = 0.9
  p = nrow(B[[1]][[1]])
  
  M = length(w)     # No. of MCMC draws
  index1 = burnin + 1
  index2 = M
  w0 = w
  w=list()
  B0 = B
  B=list()
  
  j=1
  for ( i in index1:index2){
    w[[j]] = w0[[i]]
    B[[j]] = B0[[i]]
    j = j+1
  }
  
  M = length(w)
  T = length(w[[1]])    # No. of time points
  K = length(w[[1]][[1]])     # No. of factors
  
  KM = K*M
  
  f = list()
  phi = matrix(0,KM,KM)
  L = rep(0,KM)
  
  for (m in 1:M){
      for (j in 1:K){
        l = (m-1)*K + j
        f[[l]] = rep(0,T)
        for (t in 1:T){
          f[[l]][t] = w[[m]][[t]][j]
        }
      }
  }
  
  for (l in 1:KM){
    for (m in l:KM){
      phi[l,m] = as.numeric(abs(cor(f[[l]],f[[m]])) > b_corr)
    }
    L[l] = as.numeric(sum(phi[l,]) >= M*b_draws)
  }
  
  
  
  clustr = list()
  S = which(L!=0)
  l = min(S)
  m = l:KM
  clustr[[1]] = m[phi[l,l:ncol(phi)]==1]
  C = clustr[[1]]
  c = 2
  k0 = K
  while(k0 >= (b_draws*M)){
    S = setdiff(S,C)
    if (length(S) ==0){break}
    l = min(S)
    m = setdiff(l:KM, C)
    clus_new = m[phi[l,]==1]
    k0 = length(clus_new)-sum(is.na(clus_new))
    if (k0 < (b_draws*M)){break}
    clustr[[c]] = clus_new
    C = union(C,clustr[[c]])
    c=c+1
  }
  kappa = length(clustr)
  
  f_estimate = list()
  for (c in 1:kappa){
    f_estimate[[c]] = rep(0,T)
    for (i in 1:length(clustr[[c]])){
      j = clustr[[c]][i]
      j_star = min(clustr[[c]])
      f_estimate[[c]] = f_estimate[[c]] + 
        f[[j]]* sign(cor(f[[j]],f[[j_star]]))
    }
    f_estimate[[c]] = f_estimate[[c]]/length(clustr[[c]])
  }
  
  B_new = list()
  perm_index = matrix(0,M,kappa)
  for (m in 1:M){
    f_current = list()
    for (j in 1:K){
      j0 = (m-1)*K + j
      f_current[[j]] = f[[j0]]
    }
    for (c in 1:kappa){
      corr = rep(0,K)
      for (j in 1:K){
        corr[j] = abs(cor(f_estimate[[c]],f_current[[j]]))
      }
      perm_index[m,c] = which.max(corr)
    }
    B_new[[m]] = list()
    for (t in 1:T){
      B_new[[m]][[t]] = matrix(0,p,kappa)
      for (c in 1:kappa){
        B_new[[m]][[t]][,c] = B[[m]][[t]][,perm_index[m,c]]
      }
    }
  }
  
  B_identified = Parse(B_new, burnin=1)
  B_identified
}