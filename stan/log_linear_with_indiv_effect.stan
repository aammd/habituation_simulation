data{
  int n;
  int n_tamia;
  vector[n] num_obs;
  vector[n] FID;
  array[n] int<lower=1, upper=n_tamia> tamia_id;
}
transformed data{
  vector[n] FID_log = log(FID);
}
parameters{
  real starting_FID;
  real<lower=0> sigma;
  vector[n_tamia] z;
  real<lower=0> sigma_hab;
  real<lower=0> mu_hab;
}
model{
  vector[n_tamia] habit;
  habit = z*sigma_hab + mu_hab;
  FID_log ~ normal(starting_FID + habit[tamia_id] .* num_obs, sigma);
  starting_FID ~ normal(6, 1);
  mu_hab ~ normal(-1, .5);
  z ~ std_normal();
  sigma ~ exponential(2);
}
generated quantities{
  vector[n_tamia] habit;
  habit = z*sigma_hab + mu_hab;
  vector[n] ybar = starting_FID + habit[tamia_id] .* num_obs;
  vector[n] log_lik;
  vector[n] yrep;
  for(i in 1:n){
    yrep[i] = exp(normal_rng(ybar[i], sigma));
    log_lik[i] = normal_lpdf(FID_log[i] | ybar[i], sigma);
  }
}
