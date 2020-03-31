#' Preprocess the between-host data and serialise this for stan to use.
#'
#' @param csv_file is the filename of the percentages of each strain.
#' @param dump_file is the filename to write the data to.
dumper <- function(csv_file, dump_file) {
    proportion_data <- read.csv(csv_file) / 100  # Because initially percentages.
    names(proportion_data) <- c("don", "rec")
    .d_hat_data <- proportion_data$don
    .r_hat_data <- proportion_data$rec
    .num_obs <- nrow(proportion_data)
    .sigma <- 0.05
    .sigma_tau <- 0.1
    data_list <- list(
        num_obs = .num_obs,
        d_hat = .d_hat_data,
        r_hat = .r_hat_data,
        sigma_pyro = .sigma,
        sigma_tau = .sigma_tau
    )
    saveRDS(data_list, dump_file)
}


#' Use the dumper function to prepare the data for stan.
dumper(CSV_FILE, DUMP_FILE)
