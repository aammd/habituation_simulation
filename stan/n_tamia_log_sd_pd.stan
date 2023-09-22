data{
  int n;
  int n_tamia;
  vector[n] num_obs;
  array[n] int<lower=1, upper=n_tamia> tamia_id;
  vector[n] FID;
}
transformed data{
  vector[n] ln_nobs = log(num_obs);
}
parameters{
  real logitM;
  real logitp;
  vector[n_tamia] tamia_p;
  real<lower=0> sd_logitp;
  real logd;
  vector[n_tamia] tamia_d;
  real<lower=0> sd_logd;
  real<lower=0> shape;
}
model{
  vector[n] logmu;

  logmu = - 6.9
  - log_inv_logit(logitM)
  - log1m_exp(log_inv_logit(logitp + tamia_p[tamia_id])
  + ln_nobs
  - (logd + tamia_d[tamia_id])
  - log1p_exp(ln_nobs - logd));

  FID ~ gamma(shape, shape * exp (logmu));
  tamia_p ~ normal(0, sd_logitp);
  tamia_d ~ normal(0, sd_logd);
  shape ~ cauchy(0, 5);
  logitM ~ std_normal();
  logitp ~ normal(2, .5);
  logd ~ normal(5, 2);
  sd_logitp ~ exponential(1);
}
generated quantities {
  vector[n] mu;

  mu = exp(- 6.9
  - log_inv_logit(logitM)
  - log1m_exp(log_inv_logit(logitp + tamia_p[tamia_id])
  + ln_nobs
  - (logd + tamia_d[tamia_id])
  - log1p_exp(ln_nobs - logd)));

  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = gamma_lpdf(FID[n] | shape, shape * mu);
  }
}
