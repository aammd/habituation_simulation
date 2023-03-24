# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint
# hi
# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(quarto)
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
options(clustermq.scheduler = "multicore")

options("cmdstanr_write_stan_file_dir" = "cmdstanr")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(
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
  tar_target(
    model_prior_sim,
    brm(model_bf,
        data = one_simulation,
        prior = model_priors,
        backend = "cmdstanr",
        chains = 1,
        iter = 100,
        sample_prior = "only",
        file_refit = "on_change")
  ),
  tar_target(
    prior_predictive_draws,
    command = make_prior_draws_df(brms_prior_model = model_prior_sim)
  ),
  tar_target(
    prior_predictive_df,
    command = make_unnest_prior_dataframe(
      prior_predictive_draws,
      original_data = one_simulation,
      x_name = "num_obs")
  ),
  tar_target(
    model_prior_fivetamia_sim,
    brm(model_bf,
        data = make_five_tamia(),
        prior = model_priors,
        backend = "cmdstanr",
        chains = 1,
        iter = 100,
        sample_prior = "only",
        file_refit = "on_change")
  ),
  tar_target(
    fivetamia_prior_predictive_draws,
    command = make_prior_draws_df(brms_prior_model = model_prior_fivetamia_sim)
  ),
  tar_target(
    fivetamia_prior_predictive_df,
    command = make_unnest_prior_dataframe(
      fivetamia_prior_predictive_draws,
      original_data = make_five_tamia(),
      x_name = c("num_obs", "tamia_id"))
  ),
  tar_target(
    fivetamia_group,
    fivetamia_prior_predictive_df |>
      dplyr::mutate(FID = epred) |>
      dplyr::group_by(draw_id) |>
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    demodat,
    fivetamia_prior_predictive_df |>
      dplyr::mutate(FID = epred) |>
      dplyr::filter(draw_id == 31)
  ),
  tar_target(
    sampling_model,
    command = update(model_prior_fivetamia_sim,
                     newdata = demodat,
                     sample_prior = "no",
                     chains = 2, iter = 2000)
  ),
  tar_target(
    smod2,
    command = update(sampling_model, newdata = demodat)
  ),
  tar_target(
    many_models,
    update(smod2, newdata = fivetamia_group, recompile = FALSE),
    pattern = map(fivetamia_group),
    iteration = "list"
  ),











  tar_quarto(site, path = "index.qmd")
)
