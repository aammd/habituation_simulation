data{
  int n;
  int n_tamia;
  vector[n] num_obs;
  vector[n] FID;
  array[n] int<lower=1, upper=n_tamia> tamia_id;
}
transformed data{
  vector[n] FID_sqr = sqrt(FID);
}
parameters{
  real mu_starting;
  real<lower=0> sigma_starting;
  vector[n_tamia] zs;
  real<lower=0> sigma_hab;
  real<lower=0> mu_hab;
  vector[n_tamia] zh;
  real<lower=0> sigma;
}
transformed parameters {
  // each tamia starts at a different point
  vector[n_tamia] starting_FID;
  starting_FID = zs*sigma_starting + mu_starting;
  // each tamia changes at a different rate
  vector[n_tamia] habit;
  habit = zh*sigma_hab + mu_hab;
}
model{
  // starting FID priors
  mu_starting ~ normal(6, 1);
  sigma_starting ~ exponential(1);
  zs ~ std_normal();
  // habituation priors
  mu_hab ~ normal(-1, .5);
  sigma_hab ~ exponential(1);
  zh ~ std_normal();
  // log FID is normal
  FID_sqr ~ normal(starting_FID[tamia_id] + habit[tamia_id] .* num_obs, sigma);
  sigma ~ exponential(2);
}
generated quantities{
  vector[n] ybar = starting_FID[tamia_id] + habit[tamia_id] .* num_obs;
  vector[n] log_lik;
  vector[n] yrep;
  for(i in 1:n){
    yrep[i] = exp(normal_rng(ybar[i], sigma));
    log_lik[i] = normal_lpdf(FID_sqr[i] | ybar[i], sigma);
  }
}
