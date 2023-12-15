data{
  int n;
  int n_tamia;
  vector[n] num_obs;
  vector[n] FID;
  array[n] int<lower=1,upper=n_tamia> tamia_id;
}
transformed data{
  vector[n] ln_nobs = log(num_obs);
}
parameters{
  vector[n_tamia] Z_m;
  real mu_m;
  real<lower=0> sigma_m;
  vector[n_tamia] Z_p;
  real mu_p;
  real<lower=0> sigma_p;
  vector[n_tamia] Z_d;
  real mu_d;
  real<lower=0> sigma_d;

  real<lower=0> shape;
}
transformed parameters {
  vector[n_tamia] logitM = Z_m * sigma_m + mu_m;
  vector[n_tamia] logitp = Z_p * sigma_p + mu_p;
  vector[n_tamia]   logd = Z_d * sigma_d + mu_d;
}
model{
  // mean
  vector[n] logmu;
  logmu = - 6.9
  - log_inv_logit(logitM[tamia_id])
  - log1m_exp( log_inv_logit(logitp[tamia_id])
  + ln_nobs
  - logd[tamia_id]
  - log1p_exp(ln_nobs - logd[tamia_id]) );
  // likelihood
  FID ~ gamma(shape, shape * exp (logmu));
  // priors
  Z_m ~ std_normal();
  Z_p ~ std_normal();
  Z_d ~ std_normal();
  mu_m ~ normal(1, 1);
  mu_p ~ normal(3, .5);
  mu_d ~ normal(.5, .5);
  sigma_m ~ exponential(1);
  sigma_p ~ exponential(1);
  sigma_d ~ exponential(1);
  shape ~ cauchy(0, 5);
}
generated quantities {
  vector[n] mu;
  vector[n] log_lik;
  vector[n] yrep;

  mu = exp(- 6.9
  - log_inv_logit(logitM[tamia_id])
  - log1m_exp(log_inv_logit(logitp[tamia_id])
  + ln_nobs
  - (logd[tamia_id])
  - log1p_exp(ln_nobs - logd[tamia_id])));

  for (j in 1:n) {
    log_lik[j] = gamma_lpdf(FID[j] | shape, shape * mu[j]);
    yrep[j] = gamma_rng(shape, shape * mu[j]);
  }
}

