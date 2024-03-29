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
  matrix[n_tamia, 3] Z;
  cholesky_factor_corr[3] L_Omega;
  // vector[n_tamia] Z_m;
  row_vector[3] mu_mpd;
  vector<lower=0>[3] sigma_mpd;
  // vector[n_tamia] Z_p;
  // real mu_p;
  // real<lower=0> sigma_p;
  // vector[n_tamia] Z_d;
  // real mu_d;
  // real<lower=0> sigma_d;

  real<lower=0> shape;
}
transformed parameters {
  matrix[n_tamia, 3] MPD;
  MPD = Z * diag_pre_multiply(sigma_mpd, L_Omega);
  // or just copy ar1_multilevel and create vectors? what im' doing here might not even be possible
  vector[n_tamia] logitM = mu_mpd[1] + MPD[,1];
  vector[n_tamia] logitp = mu_mpd[2] + MPD[,2];
  vector[n_tamia] logd = mu_mpd[3] + MPD[,3];
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
  to_vector(Z) ~ std_normal();
  mu_mpd[1] ~ normal(1, 1);
  mu_mpd[2] ~ normal(3, .5);
  mu_mpd[3] ~ normal(.5, .5);
  sigma_mpd ~ exponential(1);
  shape ~ lognormal(2.3, .2);
  L_Omega ~ lkj_corr_cholesky(2);
}
generated quantities {
  // for stantargets, i need to output parameters that are named the same as the simulation parameters
  real mu_m = mu_mpd[1];
  real mu_p = mu_mpd[2];
  real mu_d = mu_mpd[3];
  real sigma_m = sigma_mpd[1];
  real sigma_p = sigma_mpd[2];
  real sigma_d = sigma_mpd[3];

  vector[n] mu;
  vector[n] log_lik;
  vector[n] yrep;

  mu = exp(
    - 6.9
    - log_inv_logit(logitM[tamia_id])
    - log1m_exp(
      log_inv_logit(logitp[tamia_id])
      + ln_nobs
      - logd[tamia_id]
      - log1p_exp(ln_nobs - logd[tamia_id])
      ));

      for (j in 1:n) {
        log_lik[j] = gamma_lpdf(FID[j] | shape, shape * mu[j]);
        yrep[j] = gamma_rng(shape, shape * mu[j]);
      }
}

