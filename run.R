#!/usr/bin/env Rscript

# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

targets::tar_make()
tartargets::tar_make_clustermq(all_models_design, workers = 2) # nolint
# targets::tar_make_future(workers = 2) # nolint


targets::tar_make_future(workers = 2)
