library(xlsx)
library(dplyr)


source("src/within-host/dumper_helper.R")
#' Provides functions:
#'
#'   write_dump_file
#'   write_obs_array
#'   read_sheet
#'   read_data
#'


col_names <- names(read_sheet(DATA_FILE, "TCID50", row_numbers = row_numbers_in_sheet(DATA_FILE, "TCID50")))

tcid_data <- read_data(DATA_FILE, "TCID50")
rna_data <- read_data(DATA_FILE, "RealTime")
pyro_data <- read_data(DATA_FILE, "PYRO") / 100 # because it uses percentages.


#' Dump the observations as an array for later visualisation
write_obs_array(
    tcid_data,
    rna_data,
    pyro_data,
    gsub(".dump", "obs_array.rds", DUMP_FILE)
)

#' Dump the data so that it can be used by Stan.
write_dump_file(
    tcid_data,
    rna_data,
    pyro_data,
    col_names,
    DUMP_FILE,
    relaxation = RELAX_FACT
)
