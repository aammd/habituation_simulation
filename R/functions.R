
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



