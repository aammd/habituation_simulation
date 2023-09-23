# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint
# hi
# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(quarto)
library(stantargets)
## future for parallel computing
# library(future)
# library(future.callr)
# plan(callr)


# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("brms", "tibble",
               "tidybayes", "ggplot2",
               "tarchetypes", "dplyr"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
# options(clustermq.scheduler = "multiprocess")

options("cmdstanr_write_stan_file_dir" = here::here())

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source(files = "R")
# source("other_functions.R") # Source other scripts as needed. # nolint

## apparently the simulation parameters need to be defined OUTSIDE the pipeline?
tamia_sim_df <- expand.grid(
  max_obs = 20,
  n_tamia = c(10, 20, 40),
  logitM =  2,
  sd_logitM = 1,
  logitp = 4,
  sd_logitp = 1,
  logd = .8,
  sd_logd = 1,
  shape = 10) |>
  dplyr::mutate(sim_id = paste0("n", n_tamia))


list(
  ## starting with a simple plot of one simulation to look at the shape of the function
  tar_target(
    name = data,
    command = {
      tibble::tibble(num_obs = 0:30,
                     mu = hab_curve(num_obs,
                                    M = 200,
                                    p = .6,
                                    d = 2))
    }
  ),
  tar_target(
    name = one_curve,
    command = {
      data |>
        ggplot(aes(x = num_obs, y = mu)) +
        geom_line() +
        coord_cartesian(xlim = c(0, 30), ylim = c(0, 300)) +
        theme_bw() +
        labs(x = "N observations", y = "FID") +
        geom_vline(xintercept = 2,             linewidth = 2, lty = 2, col = "blue") +
        geom_hline(yintercept = 200,           linewidth = 2, lty = 2, col = "orange") +
        geom_hline(yintercept = 200* (1 - .6), linewidth = 2, lty = 2, col = "red")
    }
  ),
  tar_target(one_simulation,
             simulate_one_tamia()),
  tar_target(one_sim_plot,
             plot_one_tamia(one_simulation)),
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
  ## Second phase: read in and process the actual design data
  tarchetypes::tar_file_read(design,
                             command = "design.csv",
                             read = readr::read_csv(!!.x)),
  tar_target(design_data,
             command = design |>
               dplyr::rename(
                 tamia_id = ID,
                 num_obs = obs_number
               ) |>
               dplyr::mutate(FID = 200)
  ),
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
  ### EXperimental!!!!!

  tar_stan_mcmc_rep_summary(
    name = single_indiv,
    stan_files = c("stan/one_tamia.stan", "stan/one_tamia_log.stan"),
    data = one_tamia_simulation(0:25, 3, 5, .8, 20), # Runs once per rep.
    batches = 10, # Number of branch targets.
    reps = 9, # Number of model reps per branch target.
    chains = 4, # Number of MCMC chains.
    parallel_chains = 4, # How many MCMC chains to run in parallel.
    iter_warmup = 2e3, # Number of MCMC warmup iterations to run.
    iter_sampling = 2e3, # Number of MCMC post-warmup iterations to run.
    summaries = list(
      # Compute posterior intervals at levels 50% and 95%.
      # The 50% intervals should cover prior predictive parameter draws
      # 50% of the time. The 95% intervals are similar.
      # We also calculate posterior medians so we can compare them
      # directly to the prior predictive draws.
      ~posterior::quantile2(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
      # We use Gelman-Rubin potential scale reduction factors to
      # assess convergence:
      rhat = ~posterior::rhat(.x)
    ),
    deployment = "worker"),

  tar_stan_mcmc_rep_summary(
    name = single_indiv_shape10,
    stan_files = c("stan/one_tamia.stan",
                   "stan/one_tamia_log.stan",
                   "stan/one_tamia_log_shape.stan"),
    data = one_tamia_simulation(0:25, 3, 5, .8, 10), # Runs once per rep.
    batches = 3, # Number of branch targets.
    reps = 2, # Number of model reps per branch target.
    chains = 4, # Number of MCMC chains.
    refresh = 2000,
    parallel_chains = 2, # How many MCMC chains to run in parallel.
    iter_warmup = 2e3, # Number of MCMC warmup iterations to run.
    iter_sampling = 2e3, # Number of MCMC post-warmup iterations to run.
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
      # We use Gelman-Rubin potential scale reduction factors to
      # assess convergence:
      rhat = ~posterior::rhat(.x)
    ),
    deployment = "worker"),
  ## does it fit well

  tar_stan_mcmc(
    name = one_sp_fit,
    stan_files = "stan/one_tamia_log.stan",
    data = one_tamia_simulation(0:25, 3, 5, .8, 10)
  ),


  ## adding variation to one parameter
  tar_stan_mcmc_rep_summary(
    name = indiv_variation,
    stan_files = c("stan/n_tamia_log_sd_p.stan"),
    data = n_tamia_simulation_sd_p(
      num_obs = 0:25, n_tamia = 42,
      logitM =  3, logitp = 4, sd_logitp = .5,
      logd = .8, shape = 10), # Runs once per rep.
    batches = 3, # Number of branch targets.
    reps = 2, # Number of model reps per branch target.
    chains = 4, # Number of MCMC chains.
    refresh = 2000,
    parallel_chains = 2, # How many MCMC chains to run in parallel.
    iter_warmup = 2e3, # Number of MCMC warmup iterations to run.
    iter_sampling = 2e3, # Number of MCMC post-warmup iterations to run.
    summaries = list(
      ~posterior::quantile2(
        .x,
        probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
        ),
      # We use Gelman-Rubin potential scale reduction factors to
      # assess convergence:
      rhat = ~posterior::rhat(.x)
    ),
    deployment = "worker"),
  ###

  ### log transformation power analysis
  tar_target(
   data_var_p,
    command = n_tamia_simulation_sd_p(
      num_obs = 0:20, n_tamia = 10,
      logitM =  3, logitp = 4, sd_logitp = .5,
      logd = .8, shape = 10)
  ),
  tar_target(
    data_indiv,
    command = n_tamia_simulation_sd_mpd(
      max_obs = 20,
      n_tamia = 30,
      logitM =  2, sd_logitM = 1,
      logitp = 4, sd_logitp = 1,
      logd = .8, sd_logd = 1,
      shape = 10, output_mu = TRUE)
  ),
  tar_target(
    plot_indiv_test,
    command = {
      with(data_indiv,
           tibble::tibble(num_obs, tamia_id, mu)) |>
        ggplot(aes(x = num_obs, y = mu, group = tamia_id)) +
        geom_line()

    }
  ),
  tar_stan_mcmc(
    name = no_indiv,
    stan_files = "stan/log_linear_no_indiv_effect.stan",
    data = data_indiv
  ),
  tar_stan_mcmc(
    name = yes_indiv,
    stan_files = "stan/log_linear_with_indiv_effect.stan",
    data = data_indiv
  ),
  tar_stan_mcmc(
    name = no_indiv_ri,
    stan_files = "stan/log_linear_no_indiv_effect_ri.stan",
    data = data_indiv
  ),
  tar_stan_mcmc(
    name = yes_indiv_ri,
    stan_files = "stan/log_linear_with_indiv_effect_ri.stan",
    data = data_indiv
  ),
  # compare the groups --
  ## model compilation
  tar_target(
    name = non_var_log,
    command = cmdstanr::cmdstan_model(stan_file = "stan/log_linear_no_indiv_effect_ri.stan")
  ),
  tar_target(
    name = oui_var_log,
    command = cmdstanr::cmdstan_model(stan_file = "stan/log_linear_with_indiv_effect_ri.stan")
  ),
  ## simulation
  tarchetypes::tar_map_rep(
    pwr_log,
    command = compare_two_models_loo(
      model1 = non_var_log,
      model2 = oui_var_log,
      names = c("non_var_log", "oui_var_log"),
      max_obs = max_obs,
      n_tamia = n_tamia,
      logitM =  logitM,
      sd_logitM = sd_logitM,
      logitp = logitp,
      sd_logitp = sd_logitp,
      logd = logd,
      sd_logd = sd_logd,
      shape = shape
      ),
    values = tamia_sim_df,
    batches = 3,
    reps = 2,
    names = tidyselect::any_of("sim_id")
  ),
  tar_target(
    fig_loo_size,
    command = plot_loo_results(pwr_log)
  ),





  tar_quarto(site, path = "index.qmd")
)
