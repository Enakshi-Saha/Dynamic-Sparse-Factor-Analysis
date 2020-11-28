##### Simulation for DFA #####
set.seed(0)
##### Set parameters #####
n_mcmc = 10000
n_burnin = 8000
n_forecast = 5
T = n_forecast + 100
p = 10

# Among k_max factors only first k many are nonzero
k_max = floor(2*p/3)
k = 3

Theta = 0.9
phi_beta = 0.99
phi_w = 0.95

delta = 0.9
n0 = rep(20,p)
d0 = rep(0.002,p)

lambda0 = 0.9
lambda1 = 10/(1-phi_beta^2)

##### Simulate Gamma0 #####

Gamma0 = matrix(0, p, k_max)
z = floor(p/k)
block_height = 2*z

for (k1 in 1:k){
  i1 = (k1-1)*z + 1
  i2 = i1 + block_height -1
  i2 = min(i2,p) * (k1 < k) + p * (k1 == k)
  Gamma0[i1:i2,k1] = 1
}



##### Simulate Gamma0 #####

Gamma = list()
T1 = floor(T/3)
T2 = T1+10
for (t in 1:T){
  Gamma[[t]] = Gamma0
}
for (t in (T1+1) : T2){
  Gamma[[t]][,3] = numeric(p)
}

##### Simulate Beta0 #####

Beta0 = matrix(2+rnorm(p*k_max,0,1), p, k_max)

##### Simulate Beta_{1:T} #####

Beta = list()

for (t  in (1:T)){
  Beta[[t]] = matrix(0,p,k_max)
  for (j in 1:p){
    
    if (t>1){
#        thetas = theta_beta(Beta[[t-1]][j,], lambda1, lambda0, phi_beta, Theta)
        Beta[[t]][j,] = phi_beta*Beta[[t-1]][j,]+ rnorm(k_max, rep(0,k_max), 0.0025*diag(k_max))
    }else{
#        thetas = theta_beta(Beta0[j,], lambda1, lambda0, phi_beta, Theta)
        Beta[[t]][j,] = phi_beta*Beta0[j,] + rnorm(k_max, rep(0,k_max), 0.0025*diag(k_max))
    }
#    Gamma[[t]][j,] = rbinom(k,1,thetas)
  }
}

for (t in 1:T){
  Beta[[t]] = Beta[[t]] * Gamma[[t]] 
}

Beta0 = Beta0 * Gamma0

##### Simulate w0, w_{1:T} #####

w0 = rnorm(k_max, 0, 1)/(1-phi_w^2)
w = list()

for (t in 1:T){
  if (t>1){
    w[[t]] = phi_w * w[[t-1]] + rnorm(k_max, 0, 1)/(1-phi_w^2)
  }else{
    w[[t]] = phi_w * w0 + rnorm(k_max, 0, 1)/(1-phi_w^2)
  }
}

##### Simulate Y: split into train_data & test_data #####

Y = matrix(0, p, T)
for (t in 1:T){
  Y[,t] = Beta[[t]] %*% w[[t]] + rnorm(p,0,1)
}

#Y = (Y-mean(as.numeric(Y)))/sd(as.numeric(Y))

train_data = Y[,1:(T-n_forecast)]
test_data = Y[,(T-n_forecast+1):T]

Beta_true = Beta
Beta0_true = Beta0

##### Plot latent factors #####

latent = list()
for (j in 1:k_max){
  latent[[j]] = vector()
  for (t in 1:T){
    latent[[j]] = append(latent[[j]], w[[t]][[j]])
  }
  plot(latent[[j]], type = "l", xlab = "Time", ylab = "latent factor")
}

##### Estimate Beta with MCMC: sample size = NITER, burn-in period = burnin #####
ptm <- proc.time()
results = MCMC_DSS(train_data, phi_w, phi_beta_init = 0.95, phiest = TRUE,
              lambda0, lambda1, Theta, delta, n0, d0, RV = FALSE, NITER = n_mcmc)
mcmc_time = proc.time() - ptm
Beta_est = Parse(results$Beta_sim, n_burnin)
w_est = Parse(results$w_sim, n_burnin)

Beta_est_processed = post_process(results$Beta_sim, results$w_sim, burnin=n_burnin)
post_process_time = proc.time() - ptm


#Beta_credible = Credible(results$Beta_sim, 500)

# In-sample Forecast
Y_est1 = Parse(results$Y_est, n_burnin)
Y_est2 = matrix(0,p,ncol(train_data))
for (t in 1:ncol(train_data)){
  Y_est2[,t] = Beta_est[[t]] %*% w_est[[t]]
}

# In-sample MSE

mse1 = mean((as.numeric(train_data-Y_est1))^2)
mse2 = mean((as.numeric(train_data-Y_est2))^2)

##### Transform to Lower Triangular #####

##### Plot Factor Loadings #####

#Beta_est = sparsify_beta(Beta_est,0.0001)

# Draw Heatmap
# dev.off()
library(RColorBrewer)
Beta1 = list()
Beta2 = list()
Beta3 = list()

for (t in 1:ncol(train_data)){
  Beta1[[t]] = abs(Beta[[t]])
  Beta2[[t]] = abs(Beta_est[[t]])
  Beta3[[t]] = abs(Beta_est_processed[[t]])
}

library(plot.matrix)
#plot(M,col= colorRampPalette(brewer.pal(8, "Blues"))(25))

plot(Beta1[[20]], xlab = "Factors", main = "True loadings (t=20)", col= colorRampPalette(brewer.pal(8, "Blues"))(25),key=NULL)
plot(Beta2[[20]], xlab = "Factors", main = "Estimates (Lower Triangular)", col= colorRampPalette(brewer.pal(8, "Blues"))(25),key=NULL)
plot(Beta3[[20]], xlab = "Factors", main = "Estimates (Post Processed)", col= colorRampPalette(brewer.pal(8, "Blues"))(25),key=NULL)

plot(Beta1[[40]], xlab = "Factors", main = "True loadings (t=40)", col= colorRampPalette(brewer.pal(8, "Blues"))(25),key=NULL)
plot(Beta2[[40]], xlab = "Factors", main = "Estimates (Lower Triangular)", col= colorRampPalette(brewer.pal(8, "Blues"))(25),key=NULL)
plot(Beta3[[40]], xlab = "Factors", main = "Estimates (Post Processed)", col= colorRampPalette(brewer.pal(8, "Blues"))(25),key=NULL)

plot(Beta1[[80]], xlab = "Factors", main = "True loadings (t=80)", col= colorRampPalette(brewer.pal(8, "Blues"))(25),key=NULL)
plot(Beta2[[80]], xlab = "Factors", main = "Estimates (Lower Triangular)", col= colorRampPalette(brewer.pal(8, "Blues"))(25),key=NULL)
plot(Beta3[[80]], xlab = "Factors", main = "Estimates (Post Processed)", col= colorRampPalette(brewer.pal(8, "Blues"))(25),key=NULL)

##### Forecasting with MCMC_DFA #####
Y_forecast = matrix(0,p,n_forecast)
Y_forecast.var = matrix(0,p,n_forecast)
data = train_data
for (i in 1:n_forecast){
  results_temp = MCMC_DSS(data, phi_w, phi_beta_init = 0.95, phiest = TRUE,
                          lambda0, lambda1, Theta, delta, n0, d0, RV = FALSE, NITER = n_mcmc)
  
  one_step = forecast_MCMC(data, results_temp, lambda1, lambda0, delta, Theta, burnin = n_burnin, RV = FALSE)
  Y_forecast[,i] = one_step$pred_mean
  Y_forecast.var[,i] = one_step$pred_variance
  
  data = cbind(data[,-1],Y_forecast[,i])
}

mse_forecast = mean((as.numeric(test_data-Y_forecast))^2)


#heatmap(varimax(Beta2[[1]])$loadings, Colv = NA, Rowv = NA, xlab = "Factors", main = "Estimated loadings", col= colorRampPalette(brewer.pal(8, "Blues"))(25))

##### Vector Autoregressive Model #####

library(vars)
library(forecast)
data = t(train_data)
model_var1 = VAR(data, p = 1)
model_var2 = VAR(data, p = 2)
result_var1 = predict(model_var1, n.ahead = n_forecast)$fcst
result_var2 = predict(model_var2, n.ahead = n_forecast)$fcst
var_names = names(result_var1)

Y_var1 = matrix(0,p,n_forecast)
Y_var2 = Y_var1
for (i in 1:p){
  var_name = var_names[i]
  Y_var1[i,] = result_var1[[var_name]][,1]
  Y_var2[i,] = result_var2[[var_name]][,1]
}
mse_var1 = mean((as.numeric(test_data - Y_var1))^2)
mse_var2 = mean((as.numeric(test_data - Y_var2))^2)
