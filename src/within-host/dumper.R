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
#' For the values that are percentages, convert them to proportions.
pyro_data <- read_data(DATA_FILE, "PYRO")
pyro_percent_mask <- pyro_data <= 100 & pyro_data >= 0
pyro_data[pyro_percent_mask] <- pyro_data[pyro_percent_mask] / 100
#' Perform a sanity check that the pyro values are all either proportions or the
#' known error codes: -1, -2, -3 or -4.
stopifnot(class(pyro_data) == "matrix")
for (ix in 1:nrow(pyro_data)) {
    for (jx in 1:ncol(pyro_data)) {
        val <- pyro_data[ix,jx]
        good_val <- (val %in% c(-1,-2,-3,-4)) | val <= 1 | val >= 0
        if (!good_val) {
            stop("Unrecognised value in pyro data!!!")
        }
    }
}

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
