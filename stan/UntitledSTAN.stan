data{
  int n;
  vector[n] num_obs;
  vector[n] FID;
}
model{
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
  mu = 1000 * m * (1 - p * num_obs / (d + num_obs));

  FID ~ gamma(shape, shape / mu);
  shape ~ cauchy(0, 5);
  logitm ~ std_normal();
  logitp ~ normal(2, .5);
  logd ~ normal(5, 2);
}
