library(rstan)
library(dplyr)

get_dump_vars <- function() {
    c(
        "init_target_cells",
        "num_ferrets",
        "num_obs",
        "is_pure",
        "num_pure_wild",
        "num_pure_mutant",
        "num_mix",
        "max_tcid_complete",
        "tcid_complete",
        "tcid_complete_ixs",
        "max_tcid_censored",
        "tcid_censored",
        "tcid_censored_ixs",
        "tcid_ixs",
        "max_rna_complete",
        "rna_complete",
        "rna_complete_ixs",
        "rna_ixs",
        "max_pyro_complete",
        "pyro_complete",
        "pyro_complete_ixs",
        "max_pyro_left_censored",
        "pyro_left_censored",
        "pyro_left_censored_ixs",
        "max_pyro_right_censored",
        "pyro_right_censored",
        "pyro_right_censored_ixs",
        "pyro_ixs",
        "sigma_tcid",
        "sigma_rna",
        "sigma_pyro"
    )
}


#' A list of the constants needed by the Petrie model in `petrie.stan`.
petrie_constants <- list(
    "init_target_cells" = 4E8,
    "sigma_tcid" = 2,
    "sigma_rna" = 1,
    "sigma_pyro" = 0.1
)


#' Return a data frame of observations for a particular measurement from an
#' observation array with columns named after ferret identifiers. The same as
#' the output of reading the XLSX files.
#'
#' @param observation_array random observation array
#' @param obs_ix is 1 for TCID, 2 for RNA and 3 for pyrosequencing
#' @param ferret_ids a vector of ferret identifiers of the pattern:
#'   `XXXX.[0-9]{2}W[0-9]{2}V(Don|Rec).F[0-9]*`
observations_as_dataframe <- function(observation_array, obs_ix, ferret_ids) {
    result <- apply(
        X=observation_array,
        MARGIN=1,
        FUN=function(x) x[obs_ix,]
    )
    result <- as.data.frame(result)
    names(result) <- ferret_ids
    return(result)
}


#' A verbose and safe version of `saveRDS`
#'
#' @param object R object to serialise.
#' @param file the file to write to.
verbose_saveRDS <- function(object, file) {
    if (file.exists(file)) {
        cat(sprintf("There is already a file: %s\n", file))
        stop()
    }
    cat(sprintf("Writing to file: %s\n", file))
    saveRDS(object, file = file)
}
        

#' Write the simulation data to a dump file meeting the requirements of `petrie2015.stan`.
#'
#' @param tcid_data a matrix of TCID measurements
#' @param rna_data a matrix of RNA measurements
#' @param pyro_data a matrix of pyrosequencing measurements (as proportions).
#' @param col_names a vector of ferret identifiers of the pattern:
#'   `XXXX.[0-9]{2}W[0-9]{2}V(Don|Rec).F[0-9]*`
#' @param output_file is the file to write the data to
#' @param relaxation is the amount to multiply the observation sigmas by.
write_dump_file <- function(tcid_data, rna_data, pyro_data, col_names, output_file, relaxation=1) {
  num_ferrets <- dim(tcid_data)[1]
  num_days <- dim(tcid_data)[2]

  #' Return the same matrix with the some of the rows left-shifted to account
  #' for a *SINGLE* missing observation at the start of the time series.
  cycle_back_selected_rows <- function(mat) {
      rows_ixs <- (1:num_ferrets)[mat[, 1] == -1]
      new_mat <- mat
      selected_block <- mat[rows_ixs, ]
      num_cols <- ncol(mat)
      if (is.matrix(selected_block)) {
          new_mat[rows_ixs, 1:(num_cols - 1)] <- selected_block[, 2:(num_cols)]
          new_mat[rows_ixs, num_cols] <- selected_block[, 1]
      } else {
          new_mat[rows_ixs, 1:(num_cols - 1)] <- selected_block[2:(num_cols)]
          new_mat[rows_ixs, num_cols] <- selected_block[1]
      }
      return(new_mat)
  }
  
  is_wild <- grepl("100W", col_names)
  is_mix <- grepl("W[2,5,8]{1}0V", col_names)
  is_mutant <- grepl("100V", col_names)
  is_pure <- rep(NaN, num_ferrets)
  is_pure[is_wild] <- -1
  is_pure[is_mix] <- 0
  is_pure[is_mutant] <- 1
  num_pure_wild <- sum(is_pure == -1)
  num_mix <- sum(is_pure == 0)
  num_pure_mutant <- sum(is_pure == 1)
  
  init_target_cells <- 4E8
  sigma_tcid <- 2 * relaxation
  sigma_rna <- 1 * relaxation
  sigma_pyro <- 0.1 * relaxation
  #' Define the missing value code and the observation censor value.
  tcid_missing_code <- -1
  tcid_censor_value <- 0.5

  #' If the TCID50 measurement has value -1 then this indicates that no
  #' measurement was taken because the ferret was unavailable. We need to create
  #' a vector for the number of observations available for each of the ferrets.
  num_obs <- rowSums(tcid_data != tcid_missing_code)
  
  #' CAREFUL: This line mutates the TCID_DATA matrix
  tcid_data <- cycle_back_selected_rows(tcid_data)
  max_tcid_complete <- max(rowSums(tcid_data > tcid_censor_value))
  max_tcid_censored <- max(rowSums((tcid_data <= tcid_censor_value) &
    (tcid_data > tcid_missing_code)))
  tcid_complete <- matrix(-1, num_ferrets, max_tcid_complete)
  tcid_complete_ixs <- matrix(-1, num_ferrets, max_tcid_complete)
  tcid_censored <- matrix(-1, num_ferrets, max_tcid_censored)
  tcid_censored_ixs <- matrix(-1, num_ferrets, max_tcid_censored)
  #' `tcid_ixs` has a column for each ferret. The first row contains the number
  #' of complete observations for each ferret, i.e., an observation above the
  #' `tcid_censor_value`. The second row contains the number of observations
  #' below or equal to the `tcid_censor_value` but not completely missing data. 
  tcid_ixs <- matrix(-1, 2, num_ferrets)
  for (ferret in 1:num_ferrets) {
    temp_mask_complete <- tcid_data[ferret, ] > tcid_censor_value
    num_complete <- sum(temp_mask_complete)
    if (num_complete > 0) {
      tcid_complete[ferret, 1:num_complete] <- tcid_data[ferret, temp_mask_complete]
      tcid_complete_ixs[ferret, 1:num_complete] <- (1:num_days)[temp_mask_complete]
    }
    tcid_ixs[1, ferret] <- num_complete
    temp_mask_censored <- (tcid_data[ferret, ] <= tcid_censor_value) & (tcid_data[ferret, ] > tcid_missing_code)
    num_censored <- sum(temp_mask_censored)
    if (num_censored > 0) {
      tcid_censored[ferret, 1:num_censored] <- tcid_data[ferret, temp_mask_censored]
      tcid_censored_ixs[ferret, 1:num_censored] <- (1:num_days)[temp_mask_censored]
    }
    tcid_ixs[2, ferret] <- num_censored
  }
  
  rna_missing_code <- -2
  rna_censor_code <- -4  # Unclear what the censor limit it so will be ignored.
  rna_complete_min <- max(c(rna_missing_code, rna_censor_code, -1))
  # CAREFUL: This line mutates the RNA_DATA matrix
  rna_data <- cycle_back_selected_rows(rna_data)
  max_rna_complete <- max(rowSums(rna_data > rna_complete_min))
  rna_complete <- matrix(-1, num_ferrets, max_rna_complete)
  rna_complete_ixs <- matrix(-1, num_ferrets, max_rna_complete)
  rna_ixs <- rep(-1, num_ferrets)
  for (ferret in 1:num_ferrets) {
    temp_mask_complete <- rna_data[ferret, ] > rna_complete_min
    num_complete <- sum(temp_mask_complete)
    rna_complete[ferret, 1:num_complete] <- rna_data[ferret, temp_mask_complete]
    rna_complete_ixs[ferret, 1:num_complete] <- (1:num_days)[temp_mask_complete]
    rna_ixs[ferret] <- num_complete
  }
  
  pyro_missing_codes <- c(-3, -4)
  # CAREFUL: This line mutates the PYRO_DATA matrix
  pyro_data <- cycle_back_selected_rows(pyro_data)
  max_pyro_complete <- max(rowSums((pyro_data > 0) & (pyro_data < 1)))
  max_pyro_left_censored <- max(rowSums(pyro_data == 0))
  max_pyro_right_censored <- max(rowSums(pyro_data == 1))
  pyro_complete <- matrix(-1, num_ferrets, max_pyro_complete)
  pyro_complete_ixs <- matrix(-1, num_ferrets, max_pyro_complete)
  pyro_left_censored <- matrix(-1, num_ferrets, max_pyro_left_censored)
  pyro_left_censored_ixs <- matrix(-1, num_ferrets, max_pyro_left_censored)
  pyro_right_censored <- matrix(-1, num_ferrets, max_pyro_right_censored)
  pyro_right_censored_ixs <- matrix(-1, num_ferrets, max_pyro_right_censored)
  pyro_ixs <- matrix(-1, 3, num_ferrets)
  for (ferret in 1:num_ferrets) {
    temp_mask_complete <- (pyro_data[ferret, ] > 0) & (pyro_data[ferret, ] < 1)
    temp_mask_left <- pyro_data[ferret, ] == 0
    temp_mask_right <- pyro_data[ferret, ] == 1
    num_complete <- sum(temp_mask_complete)
    num_left <- sum(temp_mask_left)
    num_right <- sum(temp_mask_right)
    if (num_complete > 0) {
      pyro_complete[ferret, 1:num_complete] <- pyro_data[ferret, temp_mask_complete]
      pyro_complete_ixs[ferret, 1:num_complete] <- (1:num_days)[temp_mask_complete]
    }
    if (num_left > 0) {
      pyro_left_censored[ferret, 1:num_left] <- pyro_data[ferret, temp_mask_left]
      pyro_left_censored_ixs[ferret, 1:num_left] <- (1:num_days)[temp_mask_left]
    }
    if (num_right > 0) {
      pyro_right_censored[ferret, 1:num_right] <- pyro_data[ferret, temp_mask_right]
      pyro_right_censored_ixs[ferret, 1:num_right] <- (1:num_days)[temp_mask_right]
    }
    pyro_ixs[1, ferret] <- num_complete
    pyro_ixs[2, ferret] <- num_left
    pyro_ixs[3, ferret] <- num_right
  }
  
  
  dump_vars <- get_dump_vars()
  cat(sprintf("\nWriting the dump data to %s\n", output_file))
  stan_rdump(dump_vars, file = output_file)
}


#' Return a data frame containing the data from an XLSX sheet.
#'
#' @param data_file the XLSX file
#' @param sheet_name the name of the sheet "(TCID50|RealTime|PYRO)"
read_sheet <- function(data_file, sheet_name) {
    xlsx_sheet <- read.xlsx(
        data_file,
        sheetName = sheet_name,
        rowIndex = c(1, 3:13)
    )
    select(xlsx_sheet, matches("(Don|Rec)"))
}

#' Write the simulation data to a dump file meeting the requirements of `petrie2015.stan`.
#'
#' @param tcid_data a matrix of TCID measurements
#' @param rna_data a matrix of RNA measurements
#' @param pyro_data a matrix of pyrosequencing measurements
#' @param output_file is the file to write the data to
write_obs_array <- function(tcid_data, rna_data, pyro_data, output_file) {
    dims <- c(dim(tcid_data), 3)
    dim_perm <- c(1, 3, 2)
    obs_array <- c(tcid_data, rna_data, pyro_data) %>%
        array(dim = dims) %>%
        aperm(perm = dim_perm)
    verbose_saveRDS(obs_array, file = output_file)
}

