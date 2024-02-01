data{
  int n;
  int n_tamia;
  vector[n] num_obs;
  vector[n] FID;
  array[n] int<lower=1, upper=n_tamia> tamia_id;
}
transformed data{
  vector[n] FID_log;
  FID_log = log(FID);
}
parameters{
  real habit;
  real mu_starting;
  real<lower=0> sigma_starting;
  real<lower=0> sigma;
  vector[n_tamia] zs;
}
transformed parameters{
  // each tamia starts from a different point
  vector[n_tamia] starting_FID;
  starting_FID = zs*sigma_starting + mu_starting;
}
model{
  // average starting, sd of starting point, z-scores
  mu_starting ~ normal(6, 1);
  sigma_starting ~ exponential(1);
  zs ~ std_normal();
  // log FID is normal
  FID_log ~ normal(starting_FID[tamia_id] + habit*num_obs, sigma);
  sigma ~ exponential(2);
  habit ~ normal(-1, .5);
}
generated quantities{
  // calculate average
  vector[n] ybar;
  ybar= starting_FID[tamia_id] + habit.* num_obs;
  // log likelihood and yrep
  vector[n] log_lik;
  vector[n] yrep;
  for(i in 1:n){
    yrep[i] = exp(normal_rng(ybar[i], sigma));
    log_lik[i] = normal_lpdf(FID_log[i] | ybar[i], sigma);
  }
}
