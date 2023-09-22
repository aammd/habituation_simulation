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
  vector[n_tamia] tamia_pZ;
  real<lower=0> sd_logitp;
  real logd;
  real<lower=0> shape;
}
model{
  vector[n] logmu;
  vector[n] tamia_p = tamia_pZ*sd_logitp;

  logmu = - 6.9
  - log_inv_logit(logitM)
  - log1m_exp(log_inv_logit(logitp + tamia_p[tamia_id])
  + ln_nobs
  - logd
  - log1p_exp(ln_nobs - logd));

  FID ~ gamma(shape, shape * exp (logmu));
  tamia_pZ ~ std_normal();
  shape ~ cauchy(0, 5);
  logitM ~ std_normal();
  logitp ~ normal(2, .5);
  logd ~ normal(5, 2);
  sd_logitp ~ exponential(1);
}
