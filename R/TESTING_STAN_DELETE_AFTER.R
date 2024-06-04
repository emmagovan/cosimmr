#Testing STAN_nested
y = as.matrix(wolves_consumer[,c(1:2)])
pack = as.factor(wolves_consumer$Pack)
region = as.factor(wolves_consumer$Region)
formula = y ~ pack + (1|region) -1

a = lme4::glFormula(formula)
X_pack = a$X #packs (n_ind * n_packs)
X_region = t(as.matrix(a$reTrms$Zt)) # outside (regions) (n_ind * n_region)
wolves_sources = wolves[[3]][c(1,4,7),] #Small for now - need to accept different sources at some point
q = matrix(rep(1, 6), ncol = 2)
s_mean = wolves_sources[,c(3,5)]
s_sd = wolves_sources[,c(4,6)]
c_mean = wolves_discrimination[,c(2,4)]
c_sd = wolves_discrimination[,c(3,5)]

stan_dat = list(
  J = 2,
  N = 66,
  K = 3,
  L1 = 8,
  L2 = 3,
  y = y,
  q = q,
  s_mean = s_mean, 
  c_mean = c_mean,
  s_sd = s_sd,
  c_sd = c_sd,
  sigma_shape = c(1,1),
  sigma_rate = c(1,1),
  not_solo = 1,
  X_inner = X_pack,
  X_outer = X_region,
  hierarchical = 1
)



model = stan_model("inst/stan/STAN_nested.stan")

fit_opt <- rstan::optimizing(model,
                             data = stan_dat)

which_beta1 <- grep('beta1', names(fit_opt$par))

which_beta2 <- grep('beta2', names(fit_opt$par))

which_sigma <- grep('sigma', names(fit_opt$par))

beta1_start_opt <- structure(fit_opt$par[which_beta1],
                            dim  = c(stan_dat$K, stan_dat$L1))

beta2_start_opt <- structure(fit_opt$par[which_beta2],
                             dim  = c(stan_dat$K, stan_dat$L2))

sigma_raw_start_opt <- fit_opt$par[which_sigma][1:2]


fit_vb <- rstan::vb(
  model, data = stan_dat,
  algorithm = 'fullrank',
  pars = c('beta1', 'beta2', 'sigma'),
  init = list('beta1' = beta1_start_opt,
              'beta2' = beta2_start_opt,
              'log_sigma_raw' = sigma_raw_start_opt),
  tol_rel_obj = 0.00001, #convergence tolerance on the relative norm of the objective
  output_samples = 3600
)

#Want to simulate p here and return that?
#And return sigma in a sensible way
extracted_samples = rstan::extract(fit_vb)

#Want to extract all the betas in a sensible way first I think
beta1_ans = extracted_samples$beta1 # This is n_samples * K * n_covariates
beta2_ans = extracted_samples$beta2 
sigma_ans = extracted_samples$sigma


#CONVERT TO F
#f is x_inner * beta1 + x_outer * beta_2
#p should be n_ind * n_samples * K
f = array(NA, dim = c(66,3, 3600))

for(i in 1:66){
  for(k in 1:3){
    for(s in 1:3600){
  f[i,k,s] = X_pack[i,] %*% beta1_ans[s,k,] + X_region[i,] %*% beta2_ans[s,k,]
  }
  }
}

#Then each p is sum f / sum exp f

p_sample = array(NA, dim = c(66, 3600, 3))

for(j in 1:3600){
  for (n_obs in 1:66) {
    p_sample[n_obs,j, ] <- exp(f[n_obs,1:3, j]) / (sum((exp(f[n_obs,1:3, j]))))
  }
}


colMeans(p_sample[60,,])


#So now we have p for each individual i.e. each pack
#Should also get an average over each region??
#Set up some sort of loop dividing p up into each region
#Then colMeans I guess? to return 3600 samples






