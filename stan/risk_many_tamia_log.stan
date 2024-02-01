data{
  // sample posterior (1) or prior prediction only (0)?
  int<lower=0, upper=1>  sample_post;
  int n;
  int n_tamia;
  vector[n] num_obs;
  vector[n] FID;
  array[n] int<lower=1,upper=n_tamia> tamia_id;
  // which risk treatment is each tamia in
  array[n_tamia] int<lower=1,upper=3> risk_id;
}
transformed data{
  // log num_obs for convenience later
  vector[n] ln_nobs = log(num_obs);
}
parameters{
  // random effects for personality
  vector[n_tamia] Z_m;
  real mu_m;
  real<lower=0> sigma_m;
  vector[n_tamia] Z_p;
  real mu_p;
  real<lower=0> sigma_p;
  vector[n_tamia] Z_d;
  real mu_d;
  real<lower=0> sigma_d;
  // the effects of risk on the AVERAGE of each parameter
  // each element of the array corresponds to 3 risk treatments
  array[3] vector[3] risk_avg;
  // Shape parameter
  real<lower=0> shape;
}
transformed parameters {
  vector[n_tamia] logitM = Z_m * sigma_m + risk_avg[1][risk_id];
  vector[n_tamia] logitp = Z_p * sigma_p + risk_avg[2][risk_id];
  vector[n_tamia]   logd = Z_d * sigma_d + risk_avg[3][risk_id];
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
  if (sample_post ==1){
    FID ~ gamma(shape, shape * exp (logmu));
  }
  // priors
  Z_m ~ std_normal();
  Z_p ~ std_normal();
  Z_d ~ std_normal();
  // priors on constant terms = priors on averages
  risk_avg[1] ~ normal(1, 1);
  risk_avg[2] ~ normal(3, .5);
  risk_avg[3] ~ normal(.5, .5);
  // hierarchical variances parameters
  sigma_m ~ exponential(1);
  sigma_p ~ exponential(1);
  sigma_d ~ exponential(1);
  shape ~ lognormal(2.3, .2);
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

