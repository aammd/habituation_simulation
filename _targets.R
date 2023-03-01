# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint
# hi
# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("brms", "tibble", "tidybayes", "ggplot2", "tarchetypes"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

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
      tibble::tibble(nobs = 0:30,
                     mu = hab_curve(nobs, M = 200, p = .3, d = 2))
    }
  ),
  tar_target(
    name = one_curve,
    command = {
      data |>
        ggplot(aes(x = nobs, y = mu)) +
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
  tar_quarto(site, path = "index.qmd")
)
