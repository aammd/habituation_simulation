
# plot our habituation curve
hab_curve <- function(x, M, p, d){
  M * (1 - p*x/(d + x))
}


# simulate observations
simulate_one_tamia <- function(max_obs = 30,
                               true_alpha = 8,
                               true_m = 200,
                               true_p = .9,
                               true_d = 2){
  tibble(
    tamia_id = "un_tamia",
    num_obs = 0:max_obs,
    mu = hab_curve(num_obs,
                   M = true_m,
                   p = true_p,
                   d = true_d),
    FID = rgamma(n = length(num_obs),
                 shape = true_alpha,
                 rate = true_alpha/mu)
  )
}

plot_one_tamia <- function(df){
  df |>
    ggplot(aes(x = num_obs, y = FID)) +
    geom_point() +
    geom_line(aes(y = mu)) +
    coord_cartesian(ylim = c(0, 1000), xlim = c(0, 35)) +
    theme_bw()
}

make_many_tamia <- function(ntamia = 50){
  expand.grid(tamia_id = paste0("tamia", 1:ntamia),
              num_obs = 1:30,
              FID = 1)
}


## augment the simulated_data from model 1 with the predictor variables used by model 4:
add_indiv_vars <- function(risk_sim_list, full_df){
  pred_df <- full_df |>
    dplyr::select(Docility, Exploration, Number_captures, Sex)

  purrr::map(risk_sim_list, ~dplyr::bind_cols(pred_df, .x))
}



## updated functions for this process -- 18 May 2023

## updated functions for this process

#' Make a prior draws data frame
#'
#' @param brms_prior_model result of a call to brm with sample_prior = only
#' @param draw_vec a vector of numbers to pull out the draws you want
#' @param respname the quoted name of the response variable. This will be the
#'   name of the prior predicted observations in the output
#'
#' @return returns a tibble with one row per simulation. For each simulation
#'   ther eis a list-column, containing a data.frame with three columns:
#' * the original data (same in every simulation)
#' * the posterior epred, ie a simulation of the average. aka the deterministic part.
#' * the posterior predictions (simulated observations, suitable for refitting the model)
make_prior_draws_df <- function(brms_prior_model,
                                draw_vec = 25:44,
                                respname = "y"){
  ## rename data column

  original_data <-   brms_prior_model$data

  stopifnot(respname %in% names(original_data))

  original_data_x <- original_data[-which(names(original_data) == respname)]

  brms_prior_model |>
    as.data.frame(draw = draw_vec) |>
    dplyr::mutate(draw_id = draw_vec) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      simulated_data = list(dplyr::bind_cols(
        # extract model data
        original_data_x,
        # observations simulated from the prior
        y_new = brms::posterior_predict(
          brms_prior_model, draw_ids = draw_id) |> c(),
        # expected predicted values from the prior
        epred = brms::posterior_epred(
          brms_prior_model, draw_ids = draw_id) |> c()
      ) |>
        dplyr::rename(!!respname := y_new)
      )
    )
}

merge_prior_and_posterior_summary <- function(prior_df, posterior_df){

  prior_params_df <- prior_df |>
    dplyr::select(-lprior, -lp__, -simulated_data) |>
    tidyr::pivot_longer(cols = -draw_id, names_to = "variable")

  combined <- prior_params_df |>
    dplyr::left_join(posterior_df |>
                       tidyr::unnest(coefs),
                     by = c("draw_id", "variable"))

  return(combined)

}

calculate_post_coverage <- function(.prior_post_combined){
  .prior_post_combined |>
    dplyr::mutate(covered = value > q5 & value < q95) |>
    dplyr::group_by(variable) |>
    dplyr::summarise(prop_covered = mean(covered))
}



plot_posteriors_model4 <- function(modlist){

  model_fit_params <- modlist |>
    purrr::map_dfr(
      ~ .x |> gather_rvars(b_logitp_SexM,
                           b_logitp_scaleDocility,
                           b_logitp_scaleExploration,
                           b_logitM_scaleNumber_captures),
      .id = "sim")

  library(ggplot2)

  model_fit_params |>
    ggplot2::ggplot(ggplot2::aes(x= sim, dist = .value)) +
    ggplot2::geom_hline(yintercept = 0, col = "green", linewidth = 2)+
    tidybayes::stat_pointinterval() +
    ggplot2::facet_wrap(~.variable,
                        scales = "free") +
    ggplot2::theme(axis.text.x = element_blank())+
    ggplot2::labs(x = "Simulation ID", y = "Posterior distribution")

}

calculate_prop_model4 <- function(modlist){
  nonzero_effects <- purrr::map_dfr(modlist,
                                    ~ brms::as_draws(.x,
                                                     c("b_logitp_SexM",
                                                       "b_logitp_scaleDocility",
                                                       "b_logitp_scaleExploration",
                                                       "b_logitM_scaleNumber_captures")) |>
                                      posterior::summarise_draws(),
                                    .id = "sim") |>
    dplyr::mutate(p_nonzero = q95<0 | q5>0)


  nonzero_effects |>
    dplyr::group_by(variable) |>
    dplyr::summarise(prop_nonzero = mean(p_nonzero))

}

###


#' sets up a model for calculating
#'
#'For some reason the brms model wants to recompile by default. This
#'
#' @param one_sampled_model a model, probably the one that was sampled to generate the prior predictions
#' @param dataframe one dataframe (could be a posterior)
#'
#' @return
set_up_one_model <- function(one_sampled_model, dataframe){
  update(one_sampled_model,
         sample_prior = "no",
         newdata = dataframe,
         iter = 2000, chains = 4)
}


#' Create a dataframe of prior simulations.
#'
#' @param data_values a dataframe representing 'fake data', ie the design. It
#'   should look just like the data that the experiment will produce. it needs
#'   to have a column for the response, but that column won't be used
#' @param prior_spec a prior specification, ready to use in the model.
#' @param bf_object the output of `bf`, specifying the model.
#' @param draws_wanted the draws desired. Can be any sequence of numbers less that 1000
#'
#' @return This function runs the model simulations and returns a dataframe with
#'   one row for each simulation. It performs two steps: first, prior
#'   simulations from a brms model. Then it processes these into a dataframe with `make_prior_draws_df`
#'
simulate_from_prior <- function(data_values,
                                prior_spec,
                                bf_object,
                                draws_wanted = 25:49,
                                respname = "y"){

  brm_prior_sample_only <- brm(bf_object,
                               data = data_values,
                               prior = prior_spec,
                               sample_prior = "only",
                               chains = 2,
                               iter = 1000,
                               backend = "cmdstanr",
                               file_refit = "on_change")

  twenty_five_simulations <- make_prior_draws_df(
    brms_prior_model = brm_prior_sample_only,
    draw_vec = draws_wanted, respname = respname)

  return(twenty_five_simulations)
}

#' Calculate coverage of the posterior
#'
#' Composed function that combines three steps:
#' * summarize posterior
#' * merge with the prior predictive draws (these contain the true parameter values)
#' * calculate how many times the true value is in the interval
#'
#' @param prior_predictive_draws the dataframe of prior draws. result of `simulate_from_prior()`
#' @param all_models the list of all models fit to these draws.
#'
#' @return
#' @export
#'
#' @examples
calculate_coverage <- function(prior_predictive_draws, all_models){

  coverage_posteriors <- tibble::tibble(
    draw_id = prior_predictive_draws$draw_id,
    coefs = purrr::map(all_models, posterior::summarise_draws))

  prior_post_combined <- merge_prior_and_posterior_summary(
    prior_predictive_draws,
    coverage_posteriors)

  coverage <- calculate_post_coverage(prior_post_combined)

  return(coverage)

}
