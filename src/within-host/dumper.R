library(xlsx)
library(dplyr)


source("src/within-host/dumper_helper.R")
#' Provides functions:
#'
#'   write_dump_file
#'   write_obs_array
#'   read_sheet
#'


#' Return the data from a single sheet of the XLSX file as a transposed matrix.
#'
#' Note: this is a local function because it depends on the variable
#' `DATA_FILE` which is a variable as part of the calling environment.
#' 
#' @param sheet_name the name of the sheet "(TCID50|RealTime|PYRO)"
read_data <- function(sheet_name) {
    ftm <- function(df) {
        t(as.matrix(df))
    }
    ftm(read_sheet(DATA_FILE, sheet_name))
}


col_names <- names(read_sheet(DATA_FILE, "TCID50"))

tcid_data <- read_data("TCID50")
rna_data <- read_data("RealTime")
pyro_data <- read_data("PYRO") / 100


# Dump the observations as an array for later visualisation
write_obs_array(
    tcid_data,
    rna_data,
    pyro_data,
    gsub(".dump", "obs_array.rds", DUMP_FILE)
)

#' Only dump the data selected by the name mask. The `write_dump_file` function
#' is provided by `simulation_dumper.R`.
write_dump_file(
    tcid_data,
    rna_data,
    pyro_data,
    col_names,
    DUMP_FILE,
    relaxation = RELAX_FACT
)
