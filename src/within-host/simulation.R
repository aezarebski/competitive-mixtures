
source("src/within-host/simulation_helper.R")
#' Provides functions:
#'
#'   observe_ferret
#'   ferret
#'   run_experiment
#'   laboratory_book
#'   write_to_spreadsheet
#'

cat("Starting the simulation\n")

my_f_id <- "100W0VDon-F1"

my_params <- NULL

cat("Generating observations\n")

my_obs <- run_experiment(my_params, my_protocol, my_observer)

my_lab_book <- laboratory_book(list(list(my_f_id, my_obs)), my_protocol)

cat("Writing to XLSX\n")

write_to_spreadsheet(my_lab_book, "out/simulation/demo-simulation.xlsx")
