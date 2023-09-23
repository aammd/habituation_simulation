data{
  int n;
  vector[n] num_obs;
  vector[n] FID;
}
transformed data{
  vector[n] ln_nobs = log(num_obs);
}
parameters{
  real logitM;
  real logitp;
  real logd;
  real<lower=0> shape;
}
model{
  vector[n] logmu;
  logmu = - 6.9
  - log_inv_logit(logitM)
  - log1m_exp( log_inv_logit(logitp)
  + ln_nobs
  - logd
  - log1p_exp(ln_nobs - logd) );

  FID ~ gamma(shape, shape * exp (logmu));
  shape ~ cauchy(0, 5);
  logitM ~ std_normal();
  logitp ~ normal(2, .5);
  logd ~ normal(5, 2);
}
generated quantities {
  vector[n] mu;
  vector[n] log_lik;
  vector[n] yrep;

  mu = exp(- 6.9
  - log_inv_logit(logitM)
  - log1m_exp(log_inv_logit(logitp)
  + ln_nobs
  - (logd)
  - log1p_exp(ln_nobs - logd)));

  for (j in 1:n) {
    log_lik[j] = gamma_lpdf(FID[j] | shape, shape * mu[j]);
    yrep[j] = gamma_rng(shape, shape * mu[j]);
  }
}

