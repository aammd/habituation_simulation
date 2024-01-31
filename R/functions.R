
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


##################### functions for use with stantargets

one_tamia_simulation <-function(num_obs = 0:20,
                                logitM, logitp, logd, shape){

  m = plogis(logitM)
  p = plogis(logitp)
  d = exp(logd)
  mu = 1000 * m * (1 - p * num_obs / (d + num_obs))

  df <- tibble::tibble(num_obs = num_obs,
         FID = rgamma(length(num_obs), shape = shape, rate = shape / mu))

  list(n = length(df$num_obs),
       num_obs = df$num_obs,
       FID = df$FID,
       # special part for stantargets magic power
       .join_data = list(
         logitM = logitM,
         logitp = logitp,
         logd = logd,
         shape = shape
       )
       )
}

process_prop_log_trans <- function(sim_df_results){
  single_indiv_shape10 |>
    filter(variable != "lp__") |>
    mutate(.join_data = if_else(.name == "one_tamia_log_shape" & variable == "shape",
                                true = log(.join_data), false = .join_data)) |>
    group_by(variable, .name) |>
    summarize(in50 = sum(q25 < .join_data & .join_data < q75),
              in95 = sum(q2.5 < .join_data & .join_data < q97.5))
}

process_prop_output <- function(output_df){
  output_df |>
    filter(variable != "lp__") |>
    group_by(variable, .name) |>
    summarize(in50 = sum(q25 < .join_data & .join_data < q75),
              in95 = sum(q2.5 < .join_data & .join_data < q97.5))
}

#' Simulate tamia that vary in HOW MUCH they habituate
#'
#' all dama drop from the same point and the same rate, but reach different
#' final levels of fear.
#'
#' @param num_obs vector of observations for this animal
#' @param n_tamia number of tamia we're simulating
#' @param logitM how much they are scared at the beginning, as a proportion of 1000, on logit scale
#' @param logitp average proportion of that original value lost to habituation, logit scale
#' @param sd_logitp sd on logit scale of indiv variation in that
#' @param logd half-saturation, log scale.˜
#' @param shape shape parameter
#'
#' @return data ready for a Stan model plus the .join_data list that targets likes
n_tamia_simulation_sd_p <-function(num_obs = 0:20,
                                n_tamia,
                                logitM,
                                logitp,
                                sd_logitp,
                                logd,
                                shape){

  tamia_p <- rnorm(n_tamia, mean = 0, sd = sd_logitp)
  m = plogis(logitM)
  p_vec = plogis(logitp + tamia_p)
  d = exp(logd)

  df <- tidyr::expand_grid(
    num_obs = num_obs,
    tamia_id = 1:n_tamia) |>
    dplyr::mutate(
      p = p_vec[tamia_id],
      m = m,
      d = d,
      mu = 1000 * m * (1 - p * num_obs / (d + num_obs)),
      FID = rgamma(length(num_obs), shape = shape, rate = shape / mu)
  )

  list(n = length(df$num_obs),
       n_tamia = n_tamia,
       num_obs = df$num_obs,
       tamia_id = df$tamia_id,
       FID = df$FID,
       # special part for stantargets magic power
       .join_data = list(
         logitM = logitM,
         logitp = logitp,
         sd_logitp = sd_logitp,
         tamia_p = tamia_p,
         logd = logd,
         shape = shape
       )
  )
}

#' Simulate tamia that vary in HOW MUCH they habituate
#'
#' simulate _tamia_ that habituate, but each according to an individual value of
#' * m, the starting point
#' * p, the final ending point
#' * d, the halfway point between m and p
#'
#'
#' @param num_obs vector of observations for this animal
#' @param n_tamia number of tamia we're simulating
#' @param logitM how much they are scared at the beginning, as a proportion of 1000, on logit scale
#' @param logitp average proportion of that original value lost to habituation, logit scale
#' @param sd_logitp sd on logit scale of indiv variation in that
#' @param logd half-saturation, log scale.˜
#' @param shape shape parameter
#' @param sd_logitM  sd on logit scale of M
#' @param sd_logd sd on LOG scale of d
#'
#' @return data ready for a Stan model plus the .join_data list that targets likes
n_tamia_simulation_sd_mpd <- function(
    max_obs = 20,
    n_tamia,
    logitM,
    sd_logitM,
    logitp,
    sd_logitp,
    logd,
    sd_logd,
    shape,
    output_mu = FALSE){

  # individual effects
  tamia_M <- rnorm(n_tamia, mean = 0, sd = sd_logitM)
  tamia_p <- rnorm(n_tamia, mean = 0, sd = sd_logitp)
  tamia_d <- rnorm(n_tamia, mean = 0, sd = sd_logd)

  m_indiv = plogis(logitM + tamia_M)
  p_indiv = plogis(logitp + tamia_p)
  d_indiv = exp(logd + tamia_d)

  df <- tidyr::expand_grid(
    num_obs = 1:max_obs,
    tamia_id = 1:n_tamia) |>
    dplyr::mutate(
      p = p_indiv[tamia_id],
      m = m_indiv[tamia_id],
      d = d_indiv[tamia_id],
      mu = 1000 * m * (1 - p * num_obs / (d + num_obs)),
      FID = rgamma(length(num_obs), shape = shape, rate = shape / mu)
    )

  outlist <- list(n = length(df$num_obs),
       n_tamia = n_tamia,
       num_obs = df$num_obs,
       tamia_id = df$tamia_id,
       FID = df$FID,
       mu = df$mu,
       # special part for stantargets magic power
       .join_data = list(
         logitM = logitM,
         logitp = logitp,
         # sd_logitp = sd_logitp,
         # tamia_p = tamia_p,
         # tamia_M = tamia_M,
         # tamia_d = tamia_d,
         logd = logd,
         shape = shape
       )
  )

  if (isTRUE(output_mu)) {
    outlist$mu <- df$mu
  }

  return(outlist)
}


#' Simulate from hyperparameters
#'
#' made this to match the prior specifications in many_tamia_corr
#'
#' @param .max_obs
#' @param .n_tamia
#'
#' @return
#' @export
#'
#' @examples
n_tamia_sim_hyper <- function(.max_obs, .n_tamia){
  mu_m <- rnorm(1, 1, 1)
  sigma_m <- rexp(1, 1)
  mu_p <- rnorm(1, 3, .5)
  sigma_p <- rexp(1, 1)
  mu_d <- rnorm(1, .5, .5)
  sigma_d <- rexp(1, 1)
  shape <- rlnorm(1, 2.3, .2)

  output <- n_tamia_simulation_sd_mpd(
    max_obs = .max_obs, n_tamia = .n_tamia,
    logitM = mu_m, sd_logitM = sigma_m,
    logitp = mu_p, sd_logitp = sigma_d,
    logd   = mu_d, sd_logd   = sigma_d,
    shape =  shape)

  output$.join_data <- list(
    mu_m = mu_m,
    sigma_m = sigma_m,
    mu_p = mu_p,
    sigma_p = sigma_p,
    mu_d = mu_d,
    sigma_d = sigma_d,
    shape = shape
  )

  return(output)
}


simulate_on_design <- function(designdata){
  mu_m <- rnorm(1, 1, 1)
  sigma_m <- rexp(1, 1)
  mu_p <- rnorm(1, 3, .5)
  sigma_p <- rexp(1, 1)
  mu_d <- rnorm(1, .5, .5)
  sigma_d <- rexp(1, 1)
  shape <- rlnorm(1, 2.3, .2)

  # browser()

  design_sim_df <- designdata |>
    ## make tamia_id numeric
    mutate(tamia_id = as.numeric(as.factor(tamia_id))) |>
    group_by(tamia_id) |>
    filter(num_obs == max(num_obs)) |>
    rowwise() |>
    mutate(FID_list = list(
      n_tamia_simulation_sd_mpd(
        max_obs = num_obs,
        n_tamia = 1,
        logitM = mu_m, sd_logitM = sigma_m,
        logitp = mu_p, sd_logitp = sigma_p,
        logd = mu_d, sd_logd = sigma_d,
        shape =  shape)),
      FID_df = list(as_tibble(FID_list[c("num_obs", "FID", "tamia_id", "mu")]))
    ) |>
    select(-FID_list, -FID) |>
    rename(tamia_real_id = tamia_id, max_obs = num_obs) |>
    tidyr::unnest(FID_df)

  output <- list(
    n = length(design_sim_df$num_obs),
    n_tamia = max(design_sim_df$tamia_real_id),
    num_obs = design_sim_df$num_obs,
    tamia_id = design_sim_df$tamia_real_id,
    FID = design_sim_df$FID,
    mu = design_sim_df$mu
  )

  output$.join_data <- list(
    mu_m = mu_m,
    sigma_m = sigma_m,
    mu_p = mu_p,
    sigma_p = sigma_p,
    mu_d = mu_d,
    sigma_d = sigma_d,
    shape = shape
  )

  return(  output)
}


compare_two_models_loo <- function(model1,
                                   model2,
                                   names = c("m1", "m2"),
                                   ...){

  sim_data <-  n_tamia_simulation_sd_mpd(...)

  m1_samples <- model1$sample(data = sim_data, refresh = 0L, parallel_chains = 2)

  m1_loo <- m1_samples$loo()
  # targets::tar_load(some_groups)
  m2_samples <- model2$sample(data = sim_data, refresh=0L, parallel_chains = 2)
  m2_loo <- m2_samples$loo()
  # simulate_normal(50)
  loolist <- list(m1_loo, m2_loo) |> purrr::set_names(names)
  as.data.frame(loo::loo_compare(loolist)) |>
    tibble::rownames_to_column(var = "model")
}


# ddd <- one_tamia_simulation_sd_p(0:20, n_tamia = 12,
#                                  logitM = 3, logitp = 2,
#                                  sd_logitp = .5, logd = 2, shape = 10)
#
# ddd[c("num_obs", "tamia_id", "FID")] |>
#   as_tibble() |>
#   ggplot(aes(x = num_obs, y = FID)) + geom_point() + facet_wrap(~tamia_id)

# datalist <- n_tamia_simulation_sd_p(num_obs = 0:25, n_tamia = 45,
#                                     logitM = 3, logitp = 4, sd_logitp = 1,
#                                     logd = 3, shape = 20)
# # targets::tar_load(stan_no_indiv_var)
# # stan_no_indiv_var |> str()
# one_tamia_log <- cmdstanr::cmdstan_model("stan/one_tamia_log.stan")
# one_tamia_samples <- one_tamia_log$sample(
#   data = list(n = datalist$n,
#               num_obs = datalist$num_obs,
#               FID = datalist$FID),
#   parallel_chains = 2,
#   chains = 2,
#   refresh = 0L)
#
# one_tamia_samples$loo()
#
# n_tamia_log_sd_p <- cmdstanr::cmdstan_model("stan/n_tamia_log_sd_p.stan")
# n_tamia_log_sd_p_samples <- n_tamia_log_sd_p$sample(
#   data = list(n = datalist$n,
#               num_obs = datalist$num_obs,
#               FID = datalist$FID,
#               n_tamia = datalist$n_tamia, tamia_id = datalist$tamia_id),
#   parallel_chains = 2,
#   chains = 2,
#   refresh = 0L)

plot_loo_results <- function(loo_df) {

  big_diff <- loo_df |>
    group_by(sim_id, n_tamia, tar_batch, tar_rep) |>
    filter(elpd_diff == min(elpd_diff)) |>
    ungroup() |>
    ## REVERSE the sign when the _wrong_ model wins
    mutate(elpd_diff = if_else(model == "oui_var_log",
                               true = -elpd_diff,
                               false =elpd_diff ),
           elpd_low = elpd_diff - se_diff,
           elpd_hig = elpd_diff + se_diff)



  big_diff |>
    ungroup() |>
    as.data.frame() |>
    # filter(n_tamia == 10) |>
    ggplot(
      aes(
        x = n_tamia,
        y = elpd_diff,
        ymin = elpd_low,
        ymax = elpd_hig
      )
    )  +
    geom_pointrange(position = position_jitter(height = 0)) +
    facet_wrap(~n_tamia)

}

#' plot the loo results
#'
#'
#' @param loo_df data frame of loo results. needs columns named sim_id, n_tamia,
#'   tar_batch and tar_rep
#' @param best_model_name name of the model that "should" be the winner
#'
#' @return a ggplot
plot_loo_table <- function(loo_df,
                           best_model_name = "oui_var_log") {
  big_diff <- loo_df |>
    group_by(sim_id, n_tamia, tar_batch, tar_rep) |>
    filter(elpd_diff == min(elpd_diff)) |>
    ungroup() |>
    ## REVERSE the sign when the _wrong_ model wins
    mutate(elpd_diff = if_else(model == best_model_name,
                               true = -elpd_diff,
                               false =elpd_diff ),
           elpd_low = elpd_diff - se_diff,
           elpd_hig = elpd_diff + se_diff)

  big_diff |> # glimpse() |>
    ggplot(aes(x = n_tamia,
               y = elpd_diff,
               ymin = elpd_low, ymax = elpd_hig))  +
    geom_pointrange(position = position_jitter(height = 0, width = .5))
}
