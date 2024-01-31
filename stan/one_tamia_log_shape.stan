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
  real shape;
}
model{
  vector[n] logmu;
  logmu = - 6.9
  - log_inv_logit(logitM)
  - log1m_exp( log_inv_logit(logitp)
  + ln_nobs
  - logd
  - log1p_exp(ln_nobs - logd) );

  FID ~ gamma(exp(shape), exp (logmu + shape));
  shape ~ normal(2.3, .2);
  logitM ~ normal(3, .5);
  logitp ~ normal(5, .5);
  logd ~ normal(.8, .2);
}
