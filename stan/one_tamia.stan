data{
  int n;
  vector[n] num_obs;
  vector[n] FID;
}
parameters{
  real logitM;
  real logitp;
  real logd;
  real<lower=0> shape;
}
model{
  vector[n] mu;
  real m = inv_logit(logitM);
  real p = inv_logit(logitp);
  real d = exp(logd);
  mu = 1000 * m * (1 - p * num_obs ./ (d + num_obs));

  FID ~ gamma(shape, shape / mu);
  shape ~ lognormal(2.3, .2);
  logitM ~ normal(3, .5);
  logitp ~ normal(5, .5);
  logd ~ normal(.8, .2);
}
