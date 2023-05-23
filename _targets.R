# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint
# hi
# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(quarto)

## future for parallel computing
# library(future)
# library(future.callr)
# plan(callr)


# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("brms", "tibble",
               "tidybayes", "ggplot2",
               "tarchetypes"), # packages that your targets need to run
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


list(
  ## starting with a simple plot of one simulation to look at the shape of the function
  tar_target(
    name = data,
    command = {
      tibble::tibble(num_obs = 0:30,
                     mu = hab_curve(num_obs,
                                    M = 200,
                                    p = .3,
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
        labs(x = "N observations", y = "FID")
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

  tar_quarto(site, path = "index.qmd")
)
