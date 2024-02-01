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
library(patchwork)
## future for parallel computing
# library(future)
# library(future.callr)
# plan(callr)

# Set target options:
tar_option_set(
  packages = c(#"brms",
               "tibble",
               "tidybayes", "ggplot2", "purrr",
               "tarchetypes", "dplyr", "patchwork"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
# options(clustermq.scheduler= "multiprocess")

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
  # READ IN DATA  --------------------------------------
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


  # using Stan instead ----

  tar_stan_mcmc_rep_summary(
    name = comp_param,
    stan_files = c("stan/one_tamia.stan",
                   "stan/one_tamia_log.stan",
                    "stan/one_tamia_log_shape.stan"
                   ),
    data = one_tamia_simulation(1:25,
                                logitM = 3,
                                logitp = 5,
                                logd = .8,
                                shape =  10), # Runs once per rep.
    batches = 10, # Number of branch targets.
    reps = 2, # Number of model reps per branch target.
    chains = 2, # Number of MCMC chains.
    refresh = 0L,
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
  ## prior predictive simulation --------------------
  tar_stan_mcmc_rep_summary(
    name = prior_pred,
    stan_files = c("stan/one_tamia.stan",
                   "stan/one_tamia_log.stan",
                   "stan/one_tamia_log_shape.stan"),
    data = n_tamia_simulation_sd_mpd(
      max_obs = 25,n_tamia = 1,
      logitM = 3, sd_logitM = .5,
      logitp = 5, sd_logitp = .5,
      logd = .8, sd_logd = .2,
      shape =  10), # Runs once per rep.
    batches = 10, # Number of branch targets.
    reps = 2, # Number of model reps per branch target.
    chains = 4, # Number of MCMC chains.
    refresh = 0L,
    parallel_chains = 4, # How many MCMC chains to run in parallel.
    iter_warmup = 2e3, # Number of MCMC warmup iterations to run.
    iter_sampling = 2e3, # Number of MCMC post-warmup iterations to run.
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
      # We use Gelman-Rubin potential scale reduction factors to
      # assess convergence:
      rhat = ~posterior::rhat(.x)
    ),
    deployment = "worker"),
  ### validate heirarchical model-------------
   tar_stan_mcmc_rep_summary(
    name = cov_hier,
    stan_files = c("stan/many_tamia_log.stan", "stan/many_tamia_corr.stan"),
    data = n_tamia_sim_hyper(
      .max_obs = 25, .n_tamia = 30), # Runs once per rep.
    batches = 5, # Number of branch targets.
    reps = 1, # Number of model reps per branch target.
    chains = 4, # Number of MCMC chains.
    refresh = 0L,
    parallel_chains = 3, # How many MCMC chains to run in parallel.
    iter_warmup = 2e3, # Number of MCMC warmup iterations to run.
    iter_sampling = 2e3, # Number of MCMC post-warmup iterations to run.
    summaries = list(
      ~posterior::quantile2(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
      # We use Gelman-Rubin potential scale reduction factors to
      # assess convergence:
      rhat = ~posterior::rhat(.x)
    ),
    deployment = "worker"),

# risk models -------------------------------------------------------------
tar_target(
  design_list,
  make_risk_list(design_data)
),
tar_stan_mcmc(
  name = "prior",
  stan_files = c("stan/risk_many_tamia_log.stan",
                 "stan/risk_ordinal_many_tamia_log.stan"),
  data = list_modify(
    sample_post = 0,
    design_list
  ),
  # tiny sample number because its a prior
  chains = 1,
  iter_warmup = 100,
  iter_sampling = 100

),


### POWER ANALYSIS SECTION -----------------------------
## One fake dataset ------------------
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
    fig_indiv_test,
    command =
      with(data_indiv,
           tibble::tibble(num_obs, tamia_id, mu)) |>
        ggplot(aes(x = num_obs, y = mu, group = tamia_id)) +
        geom_line()
  ),
  tar_target(
    fig_indiv_test_points,
    command =
      with(data_indiv,
           tibble::tibble(num_obs, tamia_id, FID)) |>
      ggplot(aes(x = num_obs, y = FID, group = tamia_id)) +
      geom_point()
  ),
  ## simulate on design
  tar_target(
    design_sim,
    command = simulate_on_design(design_data)
  ),

  ## compile stan models
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
# here I'm running a test to compare a model with and without an individual slope.
# in every other way the two models are the same.
# the dataset tamua_sim_df contains the parameters for a data simulation which is different in each row.
# Everything to the right of an equals sign refers to a column in that spreadsheet
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
    batches = 10,
    reps = 2,
    names = tidyselect::any_of("sim_id")
  ),
##
  tar_target(
    fig_loo_size,
    command = plot_loo_results(pwr_log)
  ),
  ### Square Root section --------

  tar_target(
    name = non_var_sqr,
    command = cmdstanr::cmdstan_model(stan_file = "stan/sqr_linear_no_indiv_effect_ri.stan")
  ),
  tar_target(
    name = oui_var_sqr,
    command = cmdstanr::cmdstan_model(stan_file = "stan/sqr_linear_with_indiv_effect_ri.stan")
  ),
  ## simulation
  tarchetypes::tar_map_rep(
    pwr_sqr,
    command = compare_two_models_loo(
      model1 = non_var_sqr,
      model2 = oui_var_sqr,
      names = c("non_var_sqr",
                "oui_var_sqr"),
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
    batches = 10,
    reps = 2,
    names = tidyselect::any_of("sim_id")
  ),

  ## Nonlinear model section-----

  ## random effects model for individual tamia
  tar_stan_mcmc(
    name = many_sp_fit,
    stan_files = "stan/many_tamia_log.stan",
    data = data_indiv
  ),
  ## compile and compare a model with and without individual parameters
  tar_target(
    name = non_var_nonlin,
    command = cmdstanr::cmdstan_model(stan_file = "stan/many_tamia_log_sameslope.stan")
  ),
  tar_target(
    name = oui_var_nonlin,
    command = cmdstanr::cmdstan_model(stan_file = "stan/many_tamia_log.stan")
  ),
  ## simulation
  tarchetypes::tar_map_rep(
    pwr_nonlin,
    command = compare_two_models_loo(
      model1 = non_var_nonlin,
      model2 = oui_var_nonlin,
      names = c("non_var", "oui_var"),
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
    batches = 10,
    reps = 2,
    names = tidyselect::any_of("sim_id")
  ),
  ## plot results
  tar_target(
    fig_loo_loglin,
    command = plot_loo_table(pwr_log, best_model_name = "oui_var_log") +
      geom_hline(yintercept = 0, lty = 2) +
      labs(title = "log-transformed FID",
           y = "Difference to model without individual variation",
           x = "Number of individauls")
  ),
  tar_target(
    fig_loo_sqr,
    command = plot_loo_table(pwr_sqr, best_model_name = "oui_var_sqr") +
      # coord_cartesian(ylim = c(-8, 1))+
      geom_hline(yintercept = 0, lty = 2) +
      labs(title = "Square-root transformed FID",
           y = "Difference to model without individual variation",
           x = "Number of individauls")
  ),
  tar_target(
    fig_loo_nonlin,
    command = plot_loo_table(pwr_nonlin, best_model_name = "oui_var") +
      coord_cartesian(ylim = c(-300, 1))+
      geom_hline(yintercept = 0, lty = 2) +
      labs(title = "Nonlinear model of FID",
           y = "Difference to model without individual variation",
           x = "Number of individauls")
  ),
  tar_target(
    fig_panel,
    command = {fig_loo_loglin + fig_loo_sqr + fig_loo_nonlin}
  ),
  tar_target(
    fig_data_example,
    command = {fig_indiv_test + fig_indiv_test_points}
  ),


  tar_quarto(site, path = "index.qmd")
)

