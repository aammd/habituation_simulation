library(tidyverse)
library(targets)


tar_make(no_indiv_draws_log_linear_no_indiv_effect)
tar_make(no_indiv_summary_log_linear_no_indiv_effect)

tar_load(no_indiv_summary_log_linear_no_indiv_effect)
no_indiv_summary_log_linear_no_indiv_effect
tar_load(no_indiv_draws_log_linear_no_indiv_effect)
tar_load(data_indiv)


draws <- no_indiv_draws_log_linear_no_indiv_effect |>
  select(starts_with("yrep")) |>
  head(40) |>
  as.matrix()

max(data_indiv$FID)

bayesplot::ppc_dens_overlay(data_indiv$FID, yrep = draws) +
  coord_cartesian(xlim = c(0, 2000))

bayesplot::ppc_scatter_avg(data_indiv$FID, yrep = draws) +
  coord_cartesian(xlim = c(0, 2000))

## indiv var


tar_make(yes_indiv_draws_log_linear_with_indiv_effect)
tar_make(yes_indiv_summary_log_linear_with_indiv_effect)

tar_load(yes_indiv_summary_log_linear_with_indiv_effect)

yes_indiv_summary_log_linear_with_indiv_effect |>
  filter(!(variable |> str_detect("z")))
tar_load(yes_indiv_draws_log_linear_with_indiv_effect)
tar_load(data_var_p)


with_draws <- yes_indiv_draws_log_linear_with_indiv_effect |>
  select(starts_with("yrep")) |>
  head(40) |>
  as.matrix()



bayesplot::ppc_dens_overlay(data_indiv$FID, yrep = with_draws) +
  coord_cartesian(xlim = c(0, 2000))

bayesplot::ppc_scatter_avg(data_var_p$FID, yrep = with_draws) +
  coord_cartesian(xlim = c(0, 2000))


### Should I change the shape parameter tot he variance parameter just in case you guys
## like it would be easier in every way, except the connection to previous code

# make it
tar_make(data_indiv)
# plot it
tar_load(data_indiv)
with(data_indiv, tibble::tibble(num_obs, tamia_id, mu)) |>
  ggplot(aes(x = num_obs, y = mu, group = tamia_id)) + geom_line()

## i want to see these data

with(data_indiv, tibble::tibble(num_obs, tamia_id, FID)) |>
  ggplot(aes(x = num_obs, y = FID, group = tamia_id)) + geom_point()


## RI s


tar_make(yes_indiv_ri_draws_log_linear_with_indiv_effect_ri)
tar_make(yes_indiv_ri_summary_log_linear_with_indiv_effect_ri)


tar_make(no_indiv_ri_draws_log_linear_no_indiv_effect_ri)
tar_make(no_indiv_ri_summary_log_linear_no_indiv_effect_ri)

tar_load(no_indiv_ri_summary_log_linear_no_indiv_effect_ri)

tar_load(pwr_log)

library(tidyverse)
glimpse(pwr_log)
big_diff <- pwr_log |>
  group_by(sim_id, n_tamia, tar_batch, tar_rep) |>
  filter(elpd_diff == min(elpd_diff)) |>
  ungroup() |>
  ## REVERSE the sign when the _wrong_ model wins
  mutate(elpd_diff = if_else(model == "oui_var_log",
                             true = -elpd_diff,
                             false =elpd_diff ),
         elpd_low = elpd_diff - se_diff,
         elpd_hig = elpd_diff + se_diff)

big_diff |> glimpse() |>
  ggplot(aes(x = n_tamia, y = elpd_diff))  +
  geom_point()


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


### look and my predictions

tar_make(one_sp_fit_draws_one_tamia_log)
tar_load(one_sp_fit_draws_one_tamia_log)

tar_load(one_sp_fit_data)

some_yrep <- one_sp_fit_draws_one_tamia_log |>
  select(starts_with("yrep")) |>
  head(100) |>
  as.matrix()

bayesplot::ppc_dens_overlay(one_sp_fit_data$FID,
                            yrep = some_yrep)# +
coord_cartesian(xlim = c(0, 2.5))


##

tar_make(many_sp_fit_draws_many_tamia_log)
tar_load(many_sp_fit_draws_many_tamia_log)
tar_load(many_sp_fit_data)

some_yrep <- many_sp_fit_draws_many_tamia_log |>
  select(starts_with("yrep")) |>
  head(100) |>
  as.matrix()

bayesplot::ppc_dens_overlay(many_sp_fit_data$FID,
                            yrep = some_yrep)# +
coord_cartesian(xlim = c(0, 2.5))
