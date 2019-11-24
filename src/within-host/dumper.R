library(xlsx)
library(dplyr)
## library(futile.logger)
## flog.appender(appender.tee("competitive-mixtures.log"))

flog.debug("Sourcing within-host/dumper_helper.R")
source("src/within-host/dumper_helper.R")
#' Provides functions:
#'
#'   write_dump_file
#'   write_obs_array
#'   read_sheet
#'   read_data
#'

col_names <- names(read_sheet(DATA_FILE, "TCID50", row_numbers = row_numbers_in_sheet(DATA_FILE, "TCID50")))
flog.info(sprintf("There are measurements for %d ferrets", length(col_names)))
tcid_data <- read_data(DATA_FILE, "TCID50")
rna_data <- read_data(DATA_FILE, "RealTime")
pyro_data <- read_data(DATA_FILE, "PYRO") / 100 # because it uses percentages.

#' Perform a sanity check that the row names in the matrices match up.
flog.debug("Checking that the ferret identifiers match up between the data sets")
ferret_ids <- function(data_matrix) {
    sapply(strsplit(x = row.names(rna_data), split = ".", fixed = TRUE), function(x) paste(tail(x,2), collapse = "."))
}
ids_all_match <- all(c(ferret_ids(tcid_data) == ferret_ids(rna_data), ferret_ids(rna_data) == ferret_ids(pyro_data)))
if (ids_all_match) {
    flog.debug("Ferret identifiers all match")
} else {
    flog.error("Ferret identifiers do not match")
    stop("Error in creating dump files because ferret identifiers do not match")
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
