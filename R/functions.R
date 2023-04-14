
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
    coord_cartesian(ylim = c(0, 210), xlim = c(0, 35)) +
    theme_bw()
}

make_prior_draws_df <- function(brms_prior_model,
                             draw_vec = 25:44){
  brms_prior_model |>
    as.data.frame(draw = draw_vec) |>
    dplyr::mutate(draw_id = draw_vec) |>
    dplyr::rowwise() |>
    dplyr::mutate(epred = list(
      brms::posterior_epred(
        brms_prior_model, draw_ids = draw_id) |> c()),
      preds = list(
        brms::posterior_predict(
          brms_prior_model, draw_ids = draw_id) |> c()))
}

make_unnest_prior_dataframe <- function(prediction_df,
                                        original_data,
                                        x_name){
  prediction_df |>
    dplyr::mutate(num_obs = list(original_data[x_name])) |>
    tidyr::unnest(cols = c(epred, preds, num_obs))
}

make_five_tamia <- function(){
  expand.grid(tamia_id = paste0("tamia", 1:5),
              num_obs = c(0, 5, 15, 20, 25, 30),
              FID = 1)
}



## workflow -- update of this process. should this be a package?
make_prior_draws_df <- function(brms_prior_model,
                                draw_vec = 25:49){

  brms_prior_model |>
    as.data.frame(draw = draw_vec) |>
    dplyr::mutate(draw_id = draw_vec) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      simulated_data = list(dplyr::bind_cols(
        # extract model data
        brms_prior_model$data,
        # observations simulated from the prior
        posterior_predict = brms::posterior_predict(
          brms_prior_model, draw_ids = draw_id) |> c(),
        # expected predicted values from the prior
        posterior_epred = brms::posterior_epred(
          brms_prior_model, draw_ids = draw_id) |> c()
      )
      )
    )
}

simulate_from_prior <- function(data_values, prior_spec, bf_object){

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
    draw_vec = 25:49)

  return(twenty_five_simulations)
}

