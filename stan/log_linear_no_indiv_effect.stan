data{
  int n;
  vector[n] num_obs;
  vector[n] FID;
}
transformed data{
  vector[n] FID_log = log(FID);
}
parameters{
  real habit;
  real starting_FID;
  real<lower=0> sigma;
}
model{
  FID_log ~ normal(starting_FID + habit*num_obs, sigma);
  sigma ~ exponential(2);
  habit ~ normal(-1, .5);
  starting_FID ~ normal(6, 1);
}
generated quantities{
  vector[n] ybar = starting_FID + habit*num_obs;
  vector[n] log_lik;
  vector[n] yrep;
  for(i in 1:n){
    yrep[i] = exp(normal_rng(ybar[i], sigma));
    log_lik[i] = normal_lpdf(FID_log[i] | ybar[i], sigma);
  }
}
