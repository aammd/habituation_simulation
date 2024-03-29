## brms

# brms model building -----------------------------------------------------
## define the model: formula, prior. Matches what we chose for the paper.
tar_target(
  model_bf,
  bf(FID ~ inv_logit(logitM) * 1000 * (1 - inv_logit(logitp)*num_obs/(exp(logd) + num_obs)),
     logitM ~ 1 + (1 |t| tamia_id),
     logitp ~ 1 + (1 |t| tamia_id),
     logd ~ 1 + (1 |t| tamia_id),
     nl = TRUE,
     family = Gamma(link = "identity"))
),
tar_target(
  model_priors,
  c(
    prior(lkj(2), class = "cor"),
    prior(exponential(2), class = "sd", nlpar = "logd"),
    prior(exponential(4), class = "sd", nlpar = "logitM"),
    prior(exponential(2), class = "sd", nlpar = "logitp"),
    prior(normal(1.5, .5), class = "b", nlpar = "logd"),
    prior(normal(.5,.5), class = "b", nlpar = "logitM"),
    prior(normal(1, .5), class = "b", nlpar = "logitp"),
    prior(gamma(6.25, .25), class = "shape")
  )
),
# make a simulation of many observations.
tar_target(many_tamia,
           make_many_tamia(50)),
tar_target(
  prior_simulation_manytamia,
  command = simulate_from_prior(data_values = many_tamia,
                                prior_spec = model_priors,
                                bf_object = model_bf,
                                respname = "FID")
),
## refit to all of these: generate one test model, fit to all, then calculate coverage
tar_target(
  model_brm,
  brm(model_bf,
      data = prior_simulation_manytamia$simulated_data[[1]],
      prior = model_priors,
      backend = "cmdstanr",
      sample_prior = "no",
      chains = 4, cores = 4)
),
tar_target(
  df_list,
  prior_simulation_manytamia$simulated_data
),
tar_target(
  all_models,
  command = update(model_brm,
                   newdata = df_list,
                   recompile = FALSE,
                   backend = "cmdstanr",
                   chains = 4, cores = 4),
  pattern = map(df_list),
  iteration = "list"
),
tar_target(
  coverage_manytamia,
  command = calculate_coverage(prior_simulation_manytamia, all_models)
),

### prior simulation with brms on design data -------------
## simulate from the prior on the observed data
tar_target(
  prior_simulation_design,
  command = simulate_from_prior(data_values = design_data,
                                prior_spec = model_priors,
                                bf_object = model_bf,
                                draws_wanted = 50:99,
                                respname = "FID")
),
## fit models to each
tar_target(
  df_list_design,
  prior_simulation_design$simulated_data
),
## making model parallel as a test -- should change above if it works
tar_target(
  model_brm_parallel,
  brm(model_bf,
      data = df_list_design[[1]],
      prior = model_priors,
      backend = "cmdstanr",
      sample_prior = "no",
      adapt_delta = 0.95,
      chains = 4, cores = 4)
),
tar_target(
  all_models_design,
  command = update(object = model_brm_parallel,
                   newdata = df_list_design,
                   recompile = FALSE,
                   backend = "cmdstanr",
                   chains = 4, cores = 4), # specify cores here?
  pattern = map(df_list_design),
  iteration = "list"
),
tar_target(
  coverage_design,
  command = calculate_coverage(prior_simulation_design, all_models_design)
),

## validate a model with treatment effects ------------------

# make risk a factor
tar_target(
  design_data_risk,
  command = design_data |>
    dplyr::mutate(Risk = as.factor(Risk))
),
# define model 1
tar_target(
  model_1_bf,
  bf(FID ~ inv_logit(logitM) * 1000 * (1 - inv_logit(logitp)*num_obs/(exp(logd) + num_obs)),
     logitM ~ 1 + Risk + (1 |t| tamia_id),
     logitp ~ 1 + Risk + (1 |t| tamia_id),
     logd ~   1 + Risk + (1 |t| tamia_id),
     nl = TRUE,
     family = Gamma(link = "identity"))
),
tar_target(
  prior_simulation_design_risk,
  command = simulate_from_prior(data_values = design_data_risk,
                                prior_spec = model_priors,
                                bf_object = model_1_bf,
                                draws_wanted = 50:99,
                                respname = "FID")
),
## same 3 steps as above: create a model, update, calculate coverage
tar_target(
  df_list_design_risk,
  prior_simulation_design_risk$simulated_data
),
tar_target(
  model_1_brm_parallel,
  brm(model_1_bf,
      data = df_list_design_risk[[1]],
      prior = model_priors,
      backend = "cmdstanr",
      sample_prior = "no",
      adapt_delta = 0.95,
      chains = 4, cores = 4)
),
tar_target(
  all_model_1_fits,
  command = update(object = model_1_brm_parallel,
                   newdata = df_list_design_risk,
                   recompile = FALSE,
                   backend = "cmdstanr",
                   chains = 4, cores = 4), # specify cores here?
  pattern = map(df_list_design_risk),
  iteration = "list"
),
tar_target(
  coverage_model_1,
  command = calculate_coverage(prior_simulation_design_risk, all_model_1_fits)
),
## define model 4, set it up, and fit it to all these simulated datasets as well.
# define model 1
tar_target(
  model_4_bf,
  bf(FID ~ inv_logit(logitM) * 1000 * (1 - inv_logit(logitp)*num_obs/(exp(logd) + num_obs)),
     logitM ~ 1 + Risk + scale(Number_captures) + (1 |t| tamia_id),
     logitp ~ 1 + Risk + Sex + scale(Docility) + scale(Exploration) + (1 |t| tamia_id),
     logd ~   1 + Risk + (1 |t| tamia_id),
     nl = TRUE,
     family = Gamma(link = "identity"))
),
## add the predictor variables to the whole model
tar_target(
  df_list_personality,
  add_indiv_vars(df_list_design_risk, design_data_risk)
),
tar_target(
  model_4_brm_parallel,
  brm(model_4_bf,
      data = df_list_personality[[1]],
      prior = model_priors,
      backend = "cmdstanr",
      sample_prior = "no",
      adapt_delta = 0.95,
      chains = 4, cores = 4)
),
tar_target(
  all_model_4_fits,
  command = update(object = model_4_brm_parallel,
                   newdata = df_list_personality,
                   recompile = FALSE,
                   backend = "cmdstanr",
                   chains = 4, cores = 4), # specify cores here?
  pattern = map(df_list_personality),
  iteration = "list"
),
tar_target(
  model4_posterior_plot,
  command = plot_posteriors_model4(modlist = all_model_4_fits)
),
tar_target(
  nonzero_model4_table,
  command = calculate_prop_model4(all_model_4_fits)
),
## log transformed model section --------------------
tar_target(
  form_fit_log,
  command = bf(log(FID) ~ 1 + num_obs + (1 + num_obs | tamia_id),
               family = gaussian(), center = FALSE)
),
tar_target(
  priors_log,
  command = c(prior(lkj(2),              class = "cor"),
              prior(normal(6, .5), class = "b", coef = "Intercept"),
              prior(exponential(1), class = "sd", coef = "Intercept", group = "tamia_id"),
              prior(exponential(1), class = "sd", coef = "num_obs",   group = "tamia_id"),
              prior(std_normal(),      class = "b",  coef = "num_obs")
  )
),
tar_target(
  fit_log,
  command = brm(form_fit_log,
                data        = prior_simulation_manytamia$simulated_data[[7]],
                prior       = priors_log,
                seed        = 1234,
                adapt_delta = 0.95,
                core        = 3,
                iter        = 5000,
                backend     = "cmdstanr")
),
tar_target(
  fit_log_fits,
  command = update(object = fit_log,
                   newdata = df_list_design_risk,
                   recompile = FALSE,
                   backend = "cmdstanr",
                   chains = 4, cores = 4), # specify cores here?
  pattern = map(df_list_design_risk),
  iteration = "list"
),

## fit without the effect on the slope

tar_target(
  form_fit_log_noslope,
  command = bf(log(FID) ~ 1 + num_obs + (1 | tamia_id),
               family = gaussian(), center = FALSE)
),
tar_target(
  priors_log_noslope,
  command = c(
    prior(normal(6, .5),  class = "b", coef = "Intercept"),
    prior(exponential(1), class = "sd", coef = "Intercept", group = "tamia_id"),
    prior(std_normal(),   class = "b",  coef = "num_obs")
  )
),
tar_target(
  fit_log_noslope,
  command = brm(form_fit_log_noslope,
                data        = prior_simulation_manytamia$simulated_data[[7]],
                prior       = priors_log_noslope,
                seed        = 1234,
                adapt_delta = 0.95,
                core        = 3,
                iter        = 5000,
                backend     = "cmdstanr")
),
tar_target(
  fit_log_fits_noslope,
  command = update(object = fit_log_noslope,
                   newdata = df_list_design_risk,
                   recompile = FALSE,
                   backend = "cmdstanr",
                   chains = 4, cores = 4), # specify cores here?
  pattern = map(df_list_design_risk),
  iteration = "list"
),
## add LOO to each of them
tar_target(
  log_fits_slope_loo,
  command = brms::add_criterion(fit_log_fits, "loo"),
  pattern = map(fit_log_fits),
  iteration = "list"
),
tar_target(
  log_fits_noslope_loo,
  command = brms::add_criterion(fit_log_fits_noslope, "loo"),
  pattern = map(fit_log_fits_noslope),
  iteration = "list"
),
## sqrt transformed model section ---------
tar_target(
  form_fit_sqrt,
  command = bf(sqrt(FID) ~ 1 + num_obs + (1 + num_obs | tamia_id),
               family = gaussian(), center = FALSE)
),
tar_target(
  priors_sqrt,
  command = c(
    prior(lkj(2),         class = "cor"),
    prior(exponential(1), class = "sd", coef = "Intercept", group = "tamia_id"),
    prior(exponential(1), class = "sd", coef = "num_obs",   group = "tamia_id"),
    prior(cauchy(0, 5),   class = "sigma"),
    prior(normal(22, 10), class = "b",  coef = "Intercept"),
    prior(normal(-2,  1), class = "b", coef = "num_obs")
  )
),
tar_target(
  fit_sqrt,
  command = brm(form_fit_sqrt,
                data        = prior_simulation_manytamia$simulated_data[[7]],
                prior       = priors_sqrt,
                seed        = 1234,
                adapt_delta = 0.8,
                core        = 4,
                iter        = 2000,
                backend     = "cmdstanr")
),
tar_target(
  fit_sqrt_fits,
  command = update(object = fit_sqrt,
                   newdata = df_list_design_risk,
                   recompile = FALSE,
                   backend = "cmdstanr",
                   chains = 4, cores = 4), # specify cores here?
  pattern = map(df_list_design_risk),
  iteration = "list"
),

## fit without the effect on the slope

tar_target(
  form_sqrt_noslope,
  command = bf(sqrt(FID) ~ 1 + num_obs + (1 | tamia_id),
               family = gaussian(), center = FALSE)
),
tar_target(
  priors_sqrt_noslope,
  command = c(
    prior(exponential(1), class = "sd", coef = "Intercept", group = "tamia_id"),
    prior(cauchy(0, 5),   class = "sigma"),
    prior(normal(22, 10), class = "b",  coef = "Intercept"),
    prior(normal(-2,  1), class = "b", coef = "num_obs")
  )
),
tar_target(
  fit_sqrt_noslope,
  command = brm(form_sqrt_noslope,
                data        = prior_simulation_manytamia$simulated_data[[7]],
                prior       = priors_sqrt_noslope,
                seed        = 1234,
                adapt_delta = 0.95,
                core        = 3,
                iter        = 5000,
                backend     = "cmdstanr")
),
tar_target(
  fit_sqrt_fits_noslope,
  command = update(object = fit_sqrt_noslope,
                   newdata = df_list_design_risk,
                   recompile = FALSE,
                   backend = "cmdstanr",
                   chains = 4, cores = 4), # specify cores here?
  pattern = map(df_list_design_risk),
  iteration = "list"
),
## add LOO to each of them
tar_target(
  sqrt_fits_slope_loo,
  command = brms::add_criterion(fit_sqrt_fits, "loo"),
  pattern = map(fit_sqrt_fits),
  iteration = "list"
),
tar_target(
  sqrt_fits_noslope_loo,
  command = brms::add_criterion(fit_sqrt_fits_noslope, "loo"),
  pattern = map(fit_sqrt_fits_noslope),
  iteration = "list"
),
