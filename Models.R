# Description -------------------------------------------------------------

#### ### ### ## #### ### ### ## #### ### ### ## #### ### ### ## #### ### ##

# Nonlinear & Linear Models for FID Study
# Created by Catherine Čapkun-Huot
# Last modification: 2023-07-08

#### ### ### ## #### ### ### ## #### ### ### ## #### ### ### ## #### ### ##


# Libraries ---------------------------------------------------------------

# Load packages
library(tidyverse)
library(tidybayes)
library(brms)


# Select priors -----------------------------------------------------------

# Select priors before loading data

# Create fake data to get priors
fake_id <- expand.grid(tamia_id = paste0("t", 1:40),
                       num_obs = (0:15),
                       FID = 200)

## Nonlinear model ----
# Formula of the model
nonlin_form <- bf(FID ~ inv_logit(logitM) * 1000 * (1 - inv_logit(logitp)*num_obs/(exp(logd) + num_obs)),
                  logitM ~ 1 + (1 |t| tamia_id),
                  logitp ~ 1 + (1 |t| tamia_id),
                  logd   ~ 1 + (1 |t| tamia_id),
                  nl = TRUE,
                  family = Gamma(link = "identity"))


# Get priors - we use this to know what priors we need to assign
get_prior(nonlin_form, data = fake_id)

# Assign priors
nonlin_priors <- c(prior(lkj(2), class = "cor"),
                   prior(exponential(2), class = "sd", nlpar = "logd"),
                   prior(exponential(4), class = "sd", nlpar = "logitM"),
                   prior(exponential(2), class = "sd", nlpar = "logitp"),
                   prior(normal(1.5, .5), class = "b", nlpar = "logd"),
                   prior(normal(.5,.5), class = "b", nlpar = "logitM"),
                   prior(normal(-1, .2), class = "b", nlpar = "logitp"),
                   prior(gamma(6.25, .25), class = "shape"))

# Run the model with the assigned priors by sampling ONLY from the priors (not from the data or fake data)
# The predictions from this model will tell us that the priors are OK if they roughly encompass the range of biologically possible values
nonlin_fit <- brm(nonlin_form,
                  data = fake_id,
                  prior = nonlin_priors,
                  backend = "cmdstanr",
                  file = "prior_nonlin",
                  sample_prior = "only",
                  file_refit = "on_change")

# Get output
summary(nonlin_fit)

# Plot it to see better
# Plot 12 draws of 40 chipmunks
fake_tamia_id_epred <- fake_id |>
  add_epred_draws(nonlin_fit, ndraws =12)
fake_tamia_id_epred |>
  ggplot(aes(x = num_obs, y = .epred, colour = tamia_id)) +
  geom_line() +
  facet_wrap(~.draw) +
  guides(colour = "none")

# It looks fine


## Sqrt model ----
#TBD

## Log model ----
# Formula
form_fit_log <- bf(log(FID) ~  1+ num_obs + (1 + num_obs | tamia_id),
                   family = gaussian(), center = FALSE)


get_prior(form_fit_log, data = fake_id)
# Assign priors
priors_log <- c(
  prior(normal(6.5, .5), class = "b", coef = "Intercept"),
  prior(normal(-0.5, .1), class = "b",  coef = "num_obs"),
  prior(lkj(2),              class = "cor"),
  prior(exponential(1), class = "sd", coef = "Intercept",   group = "tamia_id"),
  prior(exponential(2), class = "sd", coef = "num_obs",   group = "tamia_id")
)

job::job({
#Run model from priors only
log_fit <- brm(form_fit_log,
                  data = fake_id,
                  prior = priors_log,
                  backend = "cmdstanr",
                  file = "prior_log",
                  sample_prior = "only",
                  file_refit = "on_change")
})
# Plot it to see better
# Plot 12 draws of 40 chipmunks
fake_tamia_id_epred <- fake_id |>
  add_epred_draws(log_fit, ndraws =12)
fake_tamia_id_epred |>
  ggplot(aes(x = num_obs, y = exp(.epred), colour = tamia_id)) +
  geom_line() +
  facet_wrap(~.draw) +
  guides(colour = "none") +
  coord_cartesian(ylim = c(0, 1000))

# It looks fine

## Exp model ----

# Formula
form_fit_exp_s <- bf(FID    ~ inv_logit(logitM) * 1000 * (inv_logit(logitp) * exp(-exp(logr) * num_obs) + 1 - inv_logit(logitp)),
                     logitM ~ 1 + (1 |t| ID),
                     logitp ~ 1 + (1 |t| ID),
                     logr   ~ 1 + (1 |t| ID),
                     nl     = TRUE,
                     family = Gamma(link = "identity"))
# Assign priors
priors_nonlinear_exp_s <- c(prior(lkj(2),           class = "cor"),
                            prior(exponential(2),   class = "sd", nlpar = "logr"),
                            prior(exponential(2),   class = "sd", nlpar = "logitM"),
                            prior(exponential(2),   class = "sd", nlpar = "logitp"),
                            prior(normal(0, 0.5), class = "b",    coef = "Intercept",  nlpar = "logr"),
                            prior(normal(.5, 0.5),  class = "b",  coef = "Intercept",  nlpar = "logitM"),
                            prior(normal(1, 0.2),   class = "b",  coef = "Intercept",  nlpar = "logitp"),
                            prior(gamma(6.25, .25), class = "shape"))

# Rename id column
fake_id <- fake_id |> rename(ID= tamia_id)

# Run model 4
fit_exp_prior <- brm(form_fit_exp_s,
                     data        = fake_id ,
                     prior       = priors_nonlinear_exp_s,
                     seed        = 1234,
                     sample_prior = "only",
                     adapt_delta = 0.95,
                     core        = 3,
                     iter        = 1000,  # reduce iter for sample prior
                     backend     = "cmdstanr",
                     file        = "fit_exp_prior",
                     file_refit  = "on_change")

# Plot it to see better
# Plot 12 draws of 40 chipmunks

fake_tamia_id_epred <- fake_id |>
  add_epred_draws(fit_exp_prior, ndraws =12)
fake_tamia_id_epred |>
  ggplot(aes(x = num_obs, y = .epred, colour = ID)) +
  geom_line() +
  facet_wrap(~.draw) +
  guides(colour = "none")

# Looks fine

# Data --------------------------------------------------------------------

# Load FID data - data_FID_analyses.csv -, each row is a FID trial
data_FID <- read.csv("../Habituation_predators/Data/data_FID_analyses.csv",
                     header = TRUE,
                     sep = ",")

# Data structure
data_FID <- tibble(ID      = as.factor(data_FID$ID),
                   num_obs = as.numeric(data_FID$obs_number),
                   FID     = as.numeric(data_FID$FID + 0.01), # gamma distribution needs values > 0
                   threat  = as.factor(data_FID$Approach_speed),
                   walker  = as.factor(data_FID$Observer))

# Change threat imminence category names
levels(data_FID$threat) <- list(Low  = "20", Medium = "40", High = "60")

# Get number of captures for each focal chipmunk
# Load individual characteristics data - ind_charact.csv -, each row is an individual
data_ind <- read.csv("../Habituation_predators/Data/ind_charact.csv",
                     header = TRUE,
                     sep = ";")

# Data structure
data_ind <- tibble(ID   = as.factor(data_ind$ID),
                   capt = as.numeric(data_ind$Number_captures))

# Join data on FID with data on behaviour in one combined data set
data_all <- dplyr::left_join(data_FID, data_ind, by = "ID")


# FID Models ---------------------------------------------------------------

## Nonlinear model ----

# Formula
form_fit_nonlinear <- bf(FID    ~ inv_logit(logitM) * 1000 * (1 - inv_logit(logitp)*num_obs/(exp(logd) + num_obs)),
                         logitM ~ 1 + threat + capt + walker + (1 |t| ID),
                         logitp ~ 1 + threat + walker + (1 |t| ID),
                         logd   ~ 1 + threat + walker + (1 |t| ID),
                         nl     = TRUE,
                         family = Gamma(link = "identity"))

# Assign priors
priors_nonlinear <- c(prior(lkj(2),           class = "cor"),
                      prior(exponential(2),   class = "sd", nlpar = "logd"),
                      prior(exponential(2),   class = "sd", nlpar = "logitM"),
                      prior(exponential(2),   class = "sd", nlpar = "logitp"),
                      prior(normal(1.5, 0.5), class = "b",  coef = "Intercept",  nlpar = "logd"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logd"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logd"),
                      prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logd"),
                      prior(normal(.5, 0.5),  class = "b",  coef = "Intercept",  nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "capt",       nlpar = "logitM"),
                      prior(normal(1, 0.2),   class = "b",  coef = "Intercept",  nlpar = "logitp"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitp"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitp"),
                      prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitp"),
                      prior(gamma(6.25, .25), class = "shape"))

# Run model 1
fit_nonlinear <- brm(form_fit_nonlinear,
                     data        = data_all,
                     prior       = priors_nonlinear,
                     seed        = 1234,
                     adapt_delta = 0.95,
                     core        = 3,
                     iter        = 5000,
                     backend     = "cmdstanr",
                     file        = "fit_nonlinear",
                     file_refit  = "on_change")

# Output
summary(fit_nonlinear)

# Inspect model
#shinystan::launch_shinystan(fit_nonlinear)


## Sqrt model ----

# Formula
form_fit_sqrt <- bf(sqrt(FID) ~ 1 + num_obs + #threat + threat:num_obs + capt:num_obs + walker +
                      (1 + num_obs | ID),
                    family = gaussian(), center = FALSE)


get_prior(form_fit_sqrt, data = data_all)

# Assign priors
priors_sqrt <- c(prior(lkj(2),              class = "cor"),
                 prior(exponential(1), class = "sd", coef = "Intercept", group = "ID"),
                 prior(exponential(1), class = "sd", coef = "num_obs",   group = "ID"),
                 prior(cauchy(0, 5), class = "sigma"),
                 prior(normal(22, 10),    class = "b",  coef = "Intercept"),
                 prior(normal(-2,  1),      class = "b", coef = "num_obs")#,
                 # prior(normal(0, .1),      class = "b",  coef = "num_obs:threatHigh"),
                 # prior(normal(0, .1),      class = "b",  coef = "num_obs:threatMedium"),
                 # prior(normal(0, .1),      class = "b",  coef = "num_obs:capt"),
                 # prior(normal(0, .1),      class = "b",  coef = "walkerRHB"),
                 # prior(normal(0, .1),      class = "b",  coef = "threatHigh"),
                 # prior(normal(0, .1),      class = "b",  coef = "threatMedium")
                 )


## prior check
prior_pred_sqrt <- brm(form_fit_sqrt,
    data        = data_all,
    prior       = priors_sqrt,
    seed        = 1234,
    sample_prior = "only",
    core        = 1,
    iter        = 1000,
    backend     = "cmdstanr",
    file        = "fit_sqrt",
    file_refit  = "on_change")


prior_pred_sqrt_df <- add_epred_draws(data_all,
                                      prior_pred_sqrt,
                                      ndraws = 10)

arrange(prior_pred_sqrt_df, .draw)

prior_pred_sqrt_df |>
  filter(.draw == 1) |>
  ggplot(aes(x = num_obs, y = .epred)) + geom_point()

# Plot it
ggplot() +
  geom_line(data = prior_pred_sqrt_df |> filter(.draw == 2), aes(x = num_obs, y = (.epred)^2, group = ID)) +
  facet_wrap(~.draw)

# Plot it
ggplot() +
  geom_line(data = prior_pred_sqrt_df |> filter(.draw == 2), aes(x = num_obs, y = (.epred)^2, group = ID)) +
  facet_wrap(~.draw)





p_nonlinear

               stat_lineribbon(data = obs_pred_nonlinear, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#0066CC", .width = .95, alpha = 0.15) +
               scale_fill_brewer(palette = "Blues") +
               stat_lineribbon(data = obs_pred_nonlinear, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#004C99", .width = 0) + #lines without alpha value
               guides(colour = "none") +
               coord_cartesian(ylim = c(0,1000)) +
               facet_wrap(~ threat) +
               theme_bw() +
               theme(legend.position = "none") +
               xlab("Trial order") +
               ylab("FID (cm)") +
               ggtitle("Nonlinear model")



# Run model 2
fit_sqrt <- brm(form_fit_sqrt,
                data        = data_all,
                prior       = priors_sqrt,
                seed        = 1234,
                adapt_delta = 0.95,
                core        = 3,
                iter        = 5000,
                backend     = "cmdstanr",
                file        = "fit_sqrt",
                file_refit  = "on_change")

# Output
summary(fit_sqrt)


## Log model ----

# Formula
form_fit_log <- bf(log(FID) ~ 1 + num_obs + threat + threat:num_obs + capt:num_obs + walker + (1 + num_obs | ID),
                   family = gaussian())


# Assign priors
priors_log <- c(prior(lkj(2),              class = "cor"),
                prior(student_t(3, 0, 10), class = "sd", coef = "Intercept", group = "ID"),
                prior(student_t(3, 0, 10), class = "sd", coef = "num_obs",   group = "ID"),
                prior(normal(0, 100),      class = "b",  coef = "num_obs"),
                prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatHigh"),
                prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatMedium"),
                prior(normal(0, 0.2),      class = "b",  coef = "num_obs:capt"),
                prior(normal(0, 0.2),      class = "b",  coef = "walkerRHB"),
                prior(normal(.5, 0.5),     class = "b",  coef = "threatHigh"),
                prior(normal(0, 0.2),      class = "b",  coef = "threatMedium"))


# Run model 3
fit_log <- brm(form_fit_log,
               data        = data_all,
               prior       = priors_log,
               seed        = 1234,
               adapt_delta = 0.95,
               core        = 3,
               iter        = 5000,
               backend     = "cmdstanr",
               file        = "fit_log",
               file_refit  = "on_change")

# Output
summary(fit_log)


## Gamma model ----

# Formula
form_fit_gamma <- bf(FID ~ 1 + num_obs + threat + threat:num_obs + capt:num_obs + walker + (1 + num_obs | ID),
                     family = Gamma(link = "identity"))


# Assign priors
priors_gamma <- c(prior(lkj(2),              class = "cor"),
                  prior(student_t(3, 0, 10), class = "sd", coef = "Intercept", group = "ID"),
                  prior(student_t(3, 0, 10), class = "sd", coef = "num_obs",   group = "ID"),
                  prior(normal(0, 100),      class = "b",  coef = "num_obs"),
                  prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatHigh"),
                  prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatMedium"),
                  prior(normal(0, 0.2),      class = "b",  coef = "num_obs:capt"),
                  prior(normal(0, 0.2),      class = "b",  coef = "walkerRHB"),
                  prior(normal(.5, 0.5),     class = "b",  coef = "threatHigh"),
                  prior(normal(0, 0.2),      class = "b",  coef = "threatMedium"))


# Run model 4
fit_gamma <- brm(form_fit_gamma,
                 data        = data_all,
                 prior       = priors_gamma,
                 seed        = 1234,
                 adapt_delta = 0.95,
                 core        = 3,
                 iter        = 5000,
                 backend     = "cmdstanr",
                 file        = "fit_gamma",
                 file_refit  = "on_change")

# Output
summary(fit_gamma)


## Exponential model ----

# Formula
form_fit_exp <- bf(FID    ~ inv_logit(logitM) * 1000 * (inv_logit(logitp) * exp(-exp(logr) * num_obs) + 1 - inv_logit(logitp)),
                   logitM ~ 1 + threat + capt + walker + (1 |t| ID),
                   logitp ~ 1 + threat + walker + (1 |t| ID),
                   logr   ~ 1 + threat + walker + (1 |t| ID),
                   nl     = TRUE,
                   family = Gamma(link = "identity"))


# Assign priors
priors_nonlinear_exp <- c(prior(lkj(2),           class = "cor"),
                      prior(exponential(2),   class = "sd", nlpar = "logr"),
                      prior(exponential(2),   class = "sd", nlpar = "logitM"),
                      prior(exponential(2),   class = "sd", nlpar = "logitp"),
                      prior(normal(0, 0.5), class = "b",    coef = "Intercept",  nlpar = "logr"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logr"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logr"),
                      prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logr"),
                      prior(normal(.5, 0.5),  class = "b",  coef = "Intercept",  nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "capt",       nlpar = "logitM"),
                      prior(normal(1, 0.2),   class = "b",  coef = "Intercept",  nlpar = "logitp"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitp"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitp"),
                      prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitp"),
                      prior(gamma(6.25, .25), class = "shape"))

# Run model exp
fit_exp <- brm(form_fit_exp,
               data        = data_all ,
               prior       = priors_nonlinear_exp,
               seed        = 1234,
               adapt_delta = 0.95,
               core        = 3,
               iter        = 5000,
               backend     = "cmdstanr",
               file        = "fit_exp",
               file_refit  = "on_change")


# Test of between-individual differences ----------------------------------

## Nonlinear model ----
# (1) Without random effects for the slope parameters (p and d)
# Formula
form_fit_nonlinear_noRS <- bf(FID    ~ inv_logit(logitM) * 1000 * (1 - inv_logit(logitp)*num_obs/(exp(logd) + num_obs)),
                              logitM ~ 1 + threat + capt + walker + (1 | ID),
                              logitp ~ 1 + threat + walker,
                              logd   ~ 1 + threat + walker,
                              nl     = TRUE,
                              family = Gamma(link = "identity"))
# Assign priors
priors_nonlinear_noRS <- c(prior(exponential(2),   class = "sd", nlpar = "logitM"),
                      prior(normal(1.5, 0.5), class = "b",  coef = "Intercept",  nlpar = "logd"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logd"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logd"),
                      prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logd"),
                      prior(normal(.5, 0.5),  class = "b",  coef = "Intercept",  nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitM"),
                      prior(normal(0, 0.2),   class = "b",  coef = "capt",       nlpar = "logitM"),
                      prior(normal(1, 0.2),   class = "b",  coef = "Intercept",  nlpar = "logitp"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitp"),
                      prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitp"),
                      prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitp"),
                      prior(gamma(6.25, .25), class = "shape"))

# Run model
fit_nonlinear_noRS <- brm(form_fit_nonlinear_noRS,
                          data        = data_all,
                          prior       = priors_nonlinear_noRS,
                          seed        = 1234,
                          adapt_delta = 0.95,
                          core        = 3,
                          iter        = 5000,
                          backend     = "cmdstanr",
                          file        = "fit_nonlinear_noRS",
                          file_refit  = "on_change")

# (2) Without random effects (no random intercepts; M)
# New model without chipmunk ID as random effect for all the parameters
# Formula
form_fit_nonlinear_noRE <- bf(FID    ~ inv_logit(logitM) * 1000 * (1 - inv_logit(logitp)*num_obs/(exp(logd) + num_obs)),
                         logitM ~ 1 + threat + capt + walker,
                         logitp ~ 1 + threat + walker,
                         logd   ~ 1 + threat + walker,
                         nl     = TRUE,
                         family = Gamma(link = "identity"))

# Assign priors
priors_nonlinear_noRE <- c(prior(normal(1.5, 0.5), class = "b",  coef = "Intercept",  nlpar = "logd"),
                           prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logd"),
                           prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logd"),
                           prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logd"),
                           prior(normal(.5, 0.5),  class = "b",  coef = "Intercept",  nlpar = "logitM"),
                           prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitM"),
                           prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitM"),
                           prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitM"),
                           prior(normal(0, 0.2),   class = "b",  coef = "capt",       nlpar = "logitM"),
                           prior(normal(1, 0.2),   class = "b",  coef = "Intercept",  nlpar = "logitp"),
                           prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitp"),
                           prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitp"),
                           prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitp"),
                           prior(gamma(6.25, .25), class = "shape"))

# Run model
fit_nonlinear_noRE <- brm(form_fit_nonlinear_noRE,
                          data        = data_all,
                          prior       = priors_nonlinear_noRE,
                          seed        = 1234,
                          adapt_delta = 0.95,
                          core        = 3,
                          iter        = 5000,
                          backend     = "cmdstanr",
                          file        = "fit_nonlinear_noRE",
                          file_refit  = "on_change")

### LOO ----
loo_compare(loo(fit_nonlinear), loo(fit_nonlinear_noRS), loo(fit_nonlinear_noRE))


## Sqrt model ----

# (1) Without random slopes

# Formula
form_fit_sqrt_noRS <- bf(sqrt(FID) ~ 1 + num_obs + threat + threat:num_obs + capt:num_obs + walker + (1 | ID),
                         family = gaussian(),
                         center = FALSE)

get_prior(form_fit_nonlinear_noRS, DATA = )


# Assign priors
priors_sqrt_noRS <- c(prior(student_t(3, 0, 10), class = "sd", coef = "Intercept", group = "ID"),
                      prior(normal(0, 100),      class = "b",  coef = "num_obs"),
                      prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatHigh"),
                      prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatMedium"),
                      prior(normal(0, 0.2),      class = "b",  coef = "num_obs:capt"),
                      prior(normal(0, 0.2),      class = "b",  coef = "walkerRHB"),
                      prior(normal(.5, 0.5),     class = "b",  coef = "threatHigh"),
                      prior(normal(0, 0.2),      class = "b",  coef = "threatMedium"))


# Run model
fit_sqrt_noRS <- brm(form_fit_sqrt_noRS,
                     data        = data_all,
                     prior       = priors_sqrt_noRS,
                     seed        = 1234,
                     adapt_delta = 0.95,
                     core        = 3,
                     iter        = 5000,
                     backend     = "cmdstanr",
                     file        = "fit_sqrt_noRS",
                     file_refit  = "on_change")


# (2) Without random effects (no random intercepts)
# Formula
form_fit_sqrt_noRE <- bf(sqrt(FID) ~ 1 + num_obs + threat + threat:num_obs + capt:num_obs + walker,
                         family = gaussian())


# Assign priors
priors_sqrt_noRE <- c(prior(normal(0, 100),      class = "b",  coef = "num_obs"),
                      prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatHigh"),
                      prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatMedium"),
                      prior(normal(0, 0.2),      class = "b",  coef = "num_obs:capt"),
                      prior(normal(0, 0.2),      class = "b",  coef = "walkerRHB"),
                      prior(normal(.5, 0.5),     class = "b",  coef = "threatHigh"),
                      prior(normal(0, 0.2),      class = "b",  coef = "threatMedium"))


# Run model
fit_sqrt_noRE <- brm(form_fit_sqrt_noRE,
                     data        = data_all,
                     prior       = priors_sqrt_noRE,
                     seed        = 1234,
                     adapt_delta = 0.95,
                     core        = 3,
                     iter        = 5000,
                     backend     = "cmdstanr",
                     file        = "fit_sqrt_noRE",
                     file_refit  = "on_change")

### LOO ----
loo_compare(loo(fit_sqrt), loo(fit_sqrt_noRS), loo(fit_sqrt_noRE))


## Log model ----

# (1) Without random slopes

# Formula
form_fit_log_noRS <- bf(log(FID) ~ 1 + num_obs + threat + threat:num_obs + capt:num_obs + walker + (1 | ID),
                        family = gaussian())


# Assign priors
priors_log_noRS <- c(prior(student_t(3, 0, 10), class = "sd", coef = "Intercept", group = "ID"),
                     prior(normal(0, 100),      class = "b",  coef = "num_obs"),
                     prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatHigh"),
                     prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatMedium"),
                     prior(normal(0, 0.2),      class = "b",  coef = "num_obs:capt"),
                     prior(normal(0, 0.2),      class = "b",  coef = "walkerRHB"),
                     prior(normal(.5, 0.5),     class = "b",  coef = "threatHigh"),
                     prior(normal(0, 0.2),      class = "b",  coef = "threatMedium"))


# Run model
fit_log_noRS <- brm(form_fit_log_noRS,
                    data        = data_all,
                    prior       = priors_log_noRS,
                    seed        = 1234,
                    adapt_delta = 0.95,
                    core        = 3,
                    iter        = 5000,
                    backend     = "cmdstanr",
                    file        = "fit_log_noRS",
                    file_refit  = "on_change")

# (2) Without random effects (no random intercepts)

# Formula
form_fit_log_noRE <- bf(log(FID) ~ 1 + num_obs + threat + threat:num_obs + capt:num_obs + walker,
                        family = gaussian())


# Assign priors
priors_log_noRE <- c(prior(normal(0, 100),      class = "b",  coef = "num_obs"),
                     prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatHigh"),
                     prior(normal(0, 0.2),      class = "b",  coef = "num_obs:threatMedium"),
                     prior(normal(0, 0.2),      class = "b",  coef = "num_obs:capt"),
                     prior(normal(0, 0.2),      class = "b",  coef = "walkerRHB"),
                     prior(normal(.5, 0.5),     class = "b",  coef = "threatHigh"),
                     prior(normal(0, 0.2),      class = "b",  coef = "threatMedium"))


# Run model
fit_log_noRE <- brm(form_fit_log_noRE,
                    data        = data_all,
                    prior       = priors_log_noRE,
                    seed        = 1234,
                    adapt_delta = 0.95,
                    core        = 3,
                    iter        = 5000,
                    backend     = "cmdstanr",
                    file        = "fit_log_noRE",
                    file_refit  = "on_change")

### LOO ----

loo_compare(loo(fit_log, moment_match = TRUE), loo(fit_log_noRS, moment_match = TRUE), loo(fit_log_noRE))

# Lots of pareto_k > 7
# Explore why
loo1 <- loo(fit_log)
print(loo1)

## Exp model ----

# (1) Without random effects for the slope parameters (p and r)
# Formula
form_fit_exp_noRS <- bf(FID    ~ inv_logit(logitM) * 1000 * (inv_logit(logitp) * exp(-exp(logr) * num_obs) + 1 - inv_logit(logitp)),
                        logitM ~ 1 + threat + capt + walker + (1 | ID),
                        logitp ~ 1 + threat + walker,
                        logr   ~ 1 + threat + walker,
                        nl     = TRUE,
                        family = Gamma(link = "identity"))


# Assign priors
priors_exp_noRS <- c(prior(exponential(2),   class = "sd", nlpar = "logitM"),
                     prior(normal(0, 0.5), class = "b",    coef = "Intercept",  nlpar = "logr"),
                     prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logr"),
                     prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logr"),
                     prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logr"),
                     prior(normal(.5, 0.5),  class = "b",  coef = "Intercept",  nlpar = "logitM"),
                     prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitM"),
                     prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitM"),
                     prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitM"),
                     prior(normal(0, 0.2),   class = "b",  coef = "capt",       nlpar = "logitM"),
                     prior(normal(1, 0.2),   class = "b",  coef = "Intercept",  nlpar = "logitp"),
                     prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitp"),
                     prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitp"),
                     prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitp"),
                     prior(gamma(6.25, .25), class = "shape"))

# Run model
fit_exp_noRS <- brm(form_fit_exp_noRS,
                    data        = data_all ,
                    prior       = priors_exp_noRS,
                    seed        = 1234,
                    adapt_delta = 0.95,
                    core        = 3,
                    iter        = 5000,
                    backend     = "cmdstanr",
                    file        = "fit_exp_noRS",
                    file_refit  = "on_change")


# (2) Without random effects (no random intercepts; M)
# New model without chipmunk ID as random effect for all the parameters
# Formula
form_fit_exp_noRE <- bf(FID    ~ inv_logit(logitM) * 1000 * (inv_logit(logitp) * exp(-exp(logr) * num_obs) + 1 - inv_logit(logitp)),
                        logitM ~ 1 + threat + capt + walker,
                        logitp ~ 1 + threat + walker,
                        logr   ~ 1 + threat + walker,
                        nl     = TRUE,
                        family = Gamma(link = "identity"))

# Assign priors
priors_exp_noRE <- c(prior(normal(0, 0.5),   class = "b",    coef = "Intercept",  nlpar = "logr"),
                          prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logr"),
                          prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logr"),
                          prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logr"),
                          prior(normal(.5, 0.5),  class = "b",  coef = "Intercept",  nlpar = "logitM"),
                          prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitM"),
                          prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitM"),
                          prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitM"),
                          prior(normal(0, 0.2),   class = "b",  coef = "capt",       nlpar = "logitM"),
                          prior(normal(1, 0.2),   class = "b",  coef = "Intercept",  nlpar = "logitp"),
                          prior(normal(0, 0.2),   class = "b",  coef = "threatHigh",   nlpar = "logitp"),
                          prior(normal(0, 0.2),   class = "b",  coef = "threatMedium", nlpar = "logitp"),
                          prior(normal(0, 0.2),   class = "b",  coef = "walkerRHB",  nlpar = "logitp"),
                          prior(gamma(6.25, .25), class = "shape"))

# Run model
fit_exp_noRE <- brm(form_fit_exp_noRE,
                          data        = data_all,
                          prior       = priors_exp_noRE,
                          seed        = 1234,
                          adapt_delta = 0.95,
                          core        = 3,
                          iter        = 5000,
                          backend     = "cmdstanr",
                          file        = "fit_exp_noRE",
                          file_refit  = "on_change")

### LOO ----
loo_compare(loo(fit_exp), loo(fit_exp_noRS), loo(fit_exp_noRE), loo(fit_exp_REr))


# PPC ---------------------------------------------------------------------
library(ggpubr) # needed for ggarrange - combine plots

# Posterior predictive checks
PPC_nonlinear <- pp_check(fit_nonlinear, type = "dens_overlay", ndraws = 20) +
                 xlab("FID (cm)") +
                 ylab("Frequency") +
                 ggtitle("Nonlinear model") +
                 theme(plot.margin = margin(7, 3, 3, 20, "pt")) +
                 theme(legend.position = c(.8, .8))

PPC_exp <- pp_check(fit_exp, type = "dens_overlay", ndraws = 20) +
           xlab("FID (cm)") +
           ylab("Frequency") +
           ggtitle("Exponential model") +
           theme(plot.margin = margin(7, 3, 3, 20, "pt")) +
           theme(legend.position = c(.8, .8))

PPC_sqrt <- pp_check(fit_sqrt, type = "dens_overlay", ndraws = 20) +
                     xlab("FID (cm)") +
                     ylab("Frequency") +
                     ggtitle("Sqrt model") +
                     theme(plot.margin = margin(7, 3, 3, 20, "pt")) +
                     theme(legend.position = c(.8, .8))

PPC_log <- pp_check(fit_log, type = "dens_overlay", ndraws = 20) +
                    xlab("FID (cm)") +
                    ylab("Frequency") +
                    ggtitle("Log model") +
                    theme(plot.margin = margin(7, 3, 3, 20, "pt")) +
                    theme(legend.position = c(.8, .8))

PPC <- ggarrange(PPC_nonlinear, PPC_exp, PPC_sqrt, PPC_log, ncol = 4)



# Figure ------------------------------------------------------------------
library(modelr) # needed for data_grid


## Nonlinear model ----
# Predictions from the posterior
obs_pred_nonlinear <- add_epred_draws(data_all,
                                      fit_nonlinear,
                                      seed = 1234)

# Plot it
p_nonlinear <- ggplot() +
               geom_point(data = data_all, aes(x = num_obs, y = FID), shape = 1, size = 1) +
               stat_lineribbon(data = obs_pred_nonlinear, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#0066CC", .width = .95, alpha = 0.15) +
               scale_fill_brewer(palette = "Blues") +
               stat_lineribbon(data = obs_pred_nonlinear, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#004C99", .width = 0) + #lines without alpha value
               guides(colour = "none") +
               coord_cartesian(ylim = c(0,1000)) +
               facet_wrap(~ threat) +
               theme_bw() +
               theme(legend.position = "none") +
               xlab("Trial order") +
               ylab("FID (cm)") +
               ggtitle("Nonlinear model")

## Sqrt model ----
# Predictions from the posterior
obs_pred_sqrt <- add_epred_draws(data_all,
                                 fit_sqrt,
                                 seed = 1234)
# Plot it

# Scale of analysis
p_sqrt_anal <- ggplot() +
               stat_lineribbon(data = obs_pred_sqrt, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#0066CC", .width = .95, alpha = 0.15) +
               scale_fill_brewer(palette = "Blues") +
               stat_lineribbon(data = obs_pred_sqrt, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#004C99", .width = 0) + #lines without alpha value
               guides(colour = "none") +
               facet_wrap(~ threat) +
               theme_bw() +
               theme(legend.position = "none") +
               xlab("Trial order") +
               ylab("Sqrt(FID)") +
               ggtitle("Sqrt model - scale of analysis")

# Scale of observations
p_sqrt_obs <- ggplot() +
              geom_point(data = data_all, aes(x = num_obs, y = FID), shape = 1, size = 1) +
              stat_lineribbon(data = obs_pred_sqrt, aes(x = num_obs, y = (.epred)^2, group = ID), size = 0.5, colour = "#0066CC", .width = .95, alpha = 0.15) +
              scale_fill_brewer(palette = "Blues") +
              stat_lineribbon(data = obs_pred_sqrt, aes(x = num_obs, y = (.epred)^2, group = ID), size = 0.5, colour = "#004C99", .width = 0) + #lines without alpha value
              guides(colour = "none") +
              coord_cartesian(ylim = c(0,1000)) +
              facet_wrap(~ threat) +
              theme_bw() +
              theme(legend.position = "none") +
              xlab("Trial order") +
              ylab("FID (cm)") +
              ggtitle("Sqrt model - scale of observations")

## Log model ----
# Predictions from the posterior
obs_pred_log <- add_epred_draws(data_all,
                                fit_log,
                                seed = 1234)

# Plot it
# Scale of analysis
p_log_anal <- ggplot() +
  stat_lineribbon(data = obs_pred_log, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#0066CC", .width = .95, alpha = 0.15) +
  scale_fill_brewer(palette = "Blues") +
  stat_lineribbon(data = obs_pred_log, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#004C99", .width = 0) + #lines without alpha value
  guides(colour = "none") +
  facet_wrap(~ threat) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Trial order") +
  ylab("FID (cm)") +
  ggtitle("Log model - scale of analysis")

# Scale of observations
p_log_obs <- ggplot() +
         geom_point(data = data_all, aes(x = num_obs, y = FID), shape = 1, size = 1) +
         stat_lineribbon(data = obs_pred_log, aes(x = num_obs, y = exp(.epred), group = ID), size = 0.5, colour = "#0066CC", .width = .95, alpha = 0.15) +
         scale_fill_brewer(palette = "Blues") +
         stat_lineribbon(data = obs_pred_log, aes(x = num_obs, y = exp(.epred), group = ID), size = 0.5, colour = "#004C99", .width = 0) + #lines without alpha value
         guides(colour = "none") +
         coord_cartesian(ylim = c(0,1000)) +
         facet_wrap(~ threat) +
         theme_bw() +
         theme(legend.position = "none") +
         xlab("Trial order") +
         ylab("FID (cm)") +
         ggtitle("Log model - scale of observations")

## Exp model ----
# Predictions from the posterior
obs_pred_exp <- add_epred_draws(data_all,
                                fit_exp,
                                seed = 1234)

# Plot it
# Scale of analysis
p_exp <- ggplot() +
  geom_point(data = data_all, aes(x = num_obs, y = FID), shape = 1, size = 1) +
  stat_lineribbon(data = obs_pred_exp, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#0066CC", .width = .95, alpha = 0.15) +
  scale_fill_brewer(palette = "Blues") +
  stat_lineribbon(data = obs_pred_exp, aes(x = num_obs, y = .epred, group = ID), size = 0.5, colour = "#004C99", .width = 0) + #lines without alpha value
  guides(colour = "none") +
  coord_cartesian(ylim = c(0,1000)) +
  facet_wrap(~ threat) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Trial order") +
  ylab("FID (cm)") +
  ggtitle("Exponential model")


# All plots together
p1 <- plot_grid(p_sqrt_obs, p_sqrt_anal, p_log_obs, p_log_anal, p_nonlinear, p_exp, ncol = 2, align = "v")


# Tests -------------------------------------------------------------------

#simulate
M <- 700
p <- 0.75
d <- 2.5
curve(M * (1 - p*x/(d + x)),
      xlim = c(0, 35),
      ylim = c(0, 1000))
disp <- 200
nobs <- 15
fake_tamia <- tibble(num_obs = 1:nobs,
                     reponse_moy = M * (1 - p*num_obs/(d + num_obs)),
                     FID = rgamma(nobs, reponse_moy^2/disp^2, reponse_moy/disp^2))

fake_tamia |>
  ggplot(aes(x = num_obs, y = FID)) + geom_point()
