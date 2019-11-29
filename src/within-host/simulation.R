library(testthat)

source("src/within-host/simulation_helper.R")
#' Provides functions:
#'
#'   observe_ferret
#'   ferret
#'   run_experiment
#'   laboratory_book
#'   write_to_spreadsheet
#'

#' State a couple of expectations about the environment this code is running in.
expect_that(basename(getwd()), equals("competitive-mixtures"))
output_xlsx <- "out/simulation/demo-simulation.xlsx"
expect_true(dir.exists(dirname(output_xlsx)))
cat("Starting the simulation\n")

my_f_id <- "100W0VDon_F1"
expect_that(my_f_id, is_a("character"))

my_params <- list(
  num_target_cells = 1e8,
  v_tcid_a_0 = 1,
  v_tcid_b_0 = 3,
  rho_0 = 10,
  beta_a = 1e-5,
  k_a = 2,
  delta_a = 1,
  p_a = 1,
  c_h_a = 1,
  d_inf_a = 3.12,
  xi_a = 10,
  beta_b = 1e-5,
  k_b = 2,
  delta_b = 1,
  p_b = 1,
  c_h_b = 1,
  d_inf_b = 3.12,
  xi_b = 10
)
expect_that(my_params, is_a("list"))

my_protocol <- list(num_observations = 4)
expect_that(my_protocol, is_a("list"))

my_observer <- list(
  sigma_tcid = 100,
  sigma_rna = 100,
  sigma_pyro = 0.1
)
expect_that(my_observer, is_a("list"))


cat("Generating observations\n")


my_obs <- run_experiment(my_params, my_protocol, my_observer)
expect_that(my_obs, is_a("matrix"))
expect_that(dim(my_obs), equals(c(my_protocol$num_observations, 3)))


my_lab_book <- laboratory_book(list(list(my_f_id, my_obs)))
expect_that(my_lab_book, is_a("list"))
expect_that(my_lab_book[[1]], is_a("list"))
expect_that(my_lab_book[[1]], has_names(c("sheet", "values")))



cat("Writing to XLSX\n")

write_to_spreadsheet(my_lab_book, output_xlsx)
expect_that(file.exists(output_xlsx), is_true())
