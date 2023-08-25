library(tidyverse)

plogis(1)
plogis(-1)

simulate_one_tamia(true_m = 1000*plogis(0.06),
                   true_p = plogis(-0.05),
                   true_d = exp(1.1),
                   true_alpha = 1.85
                   ) |>
  plot_one_tamia()


# risk level is 1, 2, or 3
simulate_tamia_pars <- function(risk){

  # risk means
  # from fixedef(fit1), see manuscript repo
  risk_means <- structure(
    c(-0.11151566266108, 0.4030926408911, 0.77511291493606,
      0.06486142567499, -0.05019545731342, 1.0981644806823,
      -0.20591821491591, 0.05540001748971, 1.164680427434),
    dim = c(3L, 3L),
    dimnames = list(
      c("logitM", "logitp", "logd"),
      c("risklow", "riskMedium","riskHigh")))

  ## factors by which risk is scaled
  risk_sd_effect <- c(.8, 1, 1.2)

  ## MAP estimates of SDs
  sd_vec <- c(logitm = .27, logitp = .95, logd = .49)
  sd_vec <- sd_vec * risk_sd_effect[risk]
  ## make correlation matrix
  cor_mat <- matrix(0, 3, 3)
  cor_mat[upper.tri(cor_mat)] <- c(-0.04, 0.01, -.14)
  cor_mat <- cor_mat + t(cor_mat)
  diag(cor_mat) <- 1
  ## make vcov
  vcov_params <- diag(sd_vec) %*% cor_mat %*% diag(sd_vec)


  ## draw values for a new tamia
  pars <- MASS::mvrnorm(n = 1, mu = risk_means[,risk], Sigma = vcov_params)
  return(pars)
}


pars <- simulate_tamia_pars(risk = 3)
simulate_one_tamia(true_m = 1000*plogis(pars[1]),
                   true_p = plogis(pars[2]),
                   true_d = exp(pars[3]),
                   true_alpha = 10,
                   max_obs = 15
                   )


demo <- expand_grid(n = 1:100,
            risk = 1:3) |>
  rowwise() |>
  mutate(pars = list(simulate_tamia_pars(risk)),
         simulated_FID = list(
           simulate_one_tamia(
             true_m = 1000*plogis(pars[1]),
             true_p = plogis(pars[2]),
             true_d = exp(pars[3]),
             true_alpha = 10,
             max_obs = 15
           )
           )
         )

# process into a simple data frame
sim_df <- demo |>
  select(-pars) |>
  unnest(simulated_FID) |>
  select(-tamia_id) |>
  mutate(ID = paste0(n, "-", risk))

sim_df |>
  ggplot(aes(x = num_obs, y = mu, group = ID)) +
  geom_line() +
  facet_wrap(~risk)

sim_df |>
  ggplot(aes(x = num_obs, y = FID, group = ID)) +
  geom_point() +
  facet_wrap(~risk)



### DIRECT FROM MANUSCRIPT GITHUB
# Formula
form_fit1 <- bf(FID    ~ inv_logit(logitM) * 1000 * (1 - inv_logit(logitp)*num_obs/(exp(logd) + num_obs)),
                logitM ~ 1 + risk + (1 |t| ID),
                logitp ~ 1 + risk + (1 |t| ID),
                logd   ~ 1 + risk + (1 |t| ID),
                nl     = TRUE,
                family = Gamma(link = "identity"))

# Assign priors
priors <- c(prior(lkj(2), class = "cor"),
            prior(exponential(2),   class = "sd", nlpar = "logd"),
            prior(exponential(2),   class = "sd", nlpar = "logitM"),
            prior(exponential(2),   class = "sd", nlpar = "logitp"),
            prior(normal(1.5, .5),  class = "b",  nlpar = "logd"),
            prior(normal(.5,.5),    class = "b",  nlpar = "logitM"),
            prior(normal(-.5, .5),   class = "b",  nlpar = "logitp"),
            prior(gamma(6.25, .25), class = "shape"))

# Assign priors


# Run model 1
fit1 <- brm(form_fit1,
            data        = sim_df,
            prior       = priors,
            seed        = 1234,
            adapt_delta = 0.8,
            cores       = 4,
            iter        = 2000,
            sample_prior = "yes",
            backend     = "cmdstanr",
            file        = "fit1",
            file_refit  = "on_change")




## Second model - ID as random effect nested in risk treatments ----

form_fit2 <- bf(FID    ~ inv_logit(logitM) * 1000 * (1 - inv_logit(logitp)*num_obs/(exp(logd) + num_obs)),
                logitM ~ 1 + risk + (1 | gr(ID, by = risk, id = "t")),
                logitp ~ 1 + risk + (1 | gr(ID, by = risk, id = "t")),
                logd   ~ 1 + risk + (1 | gr(ID, by = risk, id = "t")),
                nl     = TRUE,
                family = Gamma(link = "identity"))

# Run model 2
fit2 <- brm(form_fit2,
            data        = sim_df,
            prior       = priors,
            seed        = 1234,
            adapt_delta = 0.8,
            core        = 4,
            iter        = 2000,
            backend     = "cmdstanr",
            file        = "fit2",
            file_refit  = "on_change")


fit1 <- add_criterion(fit1, "loo")
fit2 <- add_criterion(fit2, "loo")

loo_compare(fit1, fit2)

# now, c




 ## for the purposes of these simulations, we ignore correlation between the parameters and treat them as independent.


# could go back to model and remove the correlations altogether


## ANYWAY
