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
  real m = inv_logit(logitM);
  real p = inv_logit(logitp);
  real d = exp(logd);
  logmu = - 6.9
  - log_inv_logit(logitM) - log1m_exp( log_inv_logit(logitp) + ln_nobs - logd - log1p_exp(ln_nobs - logd) );

  FID ~ gamma(exp(shape), exp (logmu + shape));
  shape ~ normal(3, .5);
  logitM ~ std_normal();
  logitp ~ normal(2, .5);
  logd ~ normal(5, 2);
}
