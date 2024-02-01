library(tidyverse)
library(targets)


## ADD prior simulation figures to the document!

## Add the stan programs to the document

## Describe how the test with the males works

## R code (necessary?) for simulations to validate the ordinal models..
# not a priority

## better descriptions of the individual variation test part
# math, just the equation for the mean of each.


job::job({targets::tar_make(cov_hier)})

source("R/functions.R")

dlist <- one_tamia_simulation(1:25,
                              logitM = 3,
                              logitp = 5,
                              logd = .8,
                              shape =  10)[c("n","num_obs", "FID")]

risk_many_tamia_log <- cmdstanr::cmdstan_model("stan/risk_ordinal_many_tamia_log.stan")

tar_load(design_data)

make_risk_list <- function(dataset){
  design_tamia_num <- design_data |>
    mutate(tamia_id = as.numeric(as.factor(tamia_id)),
           risk_id = as.numeric(as.factor(Risk))) |>
    ## VERY important -- make sure they are in sequence
    arrange(tamia_id) |>
    glimpse()

  design_risk_num <- design_tamia_num |>
    select(tamia_id, risk_id) |>
    unique()

  dlist <- list(
    n = nrow(design_tamia_num),
    n_tamia = nrow(design_risk_num),
    num_obs = design_tamia_num$num_obs,
    FID = design_tamia_num$FID,
    tamia_id = design_tamia_num$tamia_id,
    risk_id = design_risk_num$risk_id)
  return(dlist)
}

risk_sample <- risk_many_tamia_log$sample(
  data = dlist, chains =1 )
library(tidybayes)

risk_sample$summary(variables = "b_risk")

risk_sample |>
  gather_rvars(b_risk[param, trt]) |>
  arrange(param)






risk_sample |>
  gather_draws(yrep[rownum], ndraws = 1) |>
  bind_cols(design_tamia_num) |>
  ggplot(aes(x = num_obs, y = .value, group = tamia_id)) +
  geom_line()



dlist
one_tamia_log_shape$sample(data =  dlist,
                           chains = 1, parallel_chains = 1)

one_tamia <- cmdstanr::cmdstan_model("stan/one_tamia.stan")
one_tamia$sample(data =  dlist,
                           chains = 1, parallel_chains = 1)



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


one_dataset <- n_tamia_sim_hyper(
      .max_obs = 25, .n_tamia = 30)

one_dataset$.join_data

many_tamia_corr <- cmdstanr::cmdstan_model("stan/many_tamia_corr.stan")
many_tamia_corr_post <- many_tamia_corr$sample(data = one_dataset, parallel_chains = 2, chains = 1)

many_tamia_corr_post$summary(variables = c("sigma_mpd", "mu_mpd", "shape"))


