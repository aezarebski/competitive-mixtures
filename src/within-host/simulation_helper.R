#' Helper functions for simulating ferret competitive mixture experiments.

#' \strong{IO}: Observe the ferret
#'
#' @param hidden_state is the internal state of the ferret from \code{ferret}
#' @param observer is the parameters of the observation model
#'
#' @return \code{IO} observations
#'
observe_ferret <- function(hidden_state, observer) {
    NULL
}


#' The hidden state of the ferret throughout an experiment
#'
#' @param params is the parameters of the ODE model
#' @param protocol is the parameters of the experimental setup
#'
#' @return hidden state of ferret
#'
ferret <- function(params, protocol) {
    NULL
}


#' \strong{IO}: Run the experiment and return the observations.
#
#' @param params is the parameters of the ODE model
#' @param protocol is the parameters of the experimental setup
#' @param observer is the parameters of the observation model
#'
#' @return \code{IO} observations
run_experiment <- function(params, protocol, observer) {
    hidden_state <- ferret(params, protocol)
    observe_experiment(hidden_state, observer)
}

#' Results of an experiment in a convenient format
#'
#' @param ferrets_and_observations lists of pairs of ferret identifiers and observations
#' @param protocol the common protocol of the experiments
#'
#' @return laboratory book
#'
laboratory_book <- function(ferrets_and_observations, protocol) {
    NULL
}

#' \strong{IO}: Write the results to a spreadsheet
#'
#' @param lab_book is a labbook from \code{laboratory_book}
#' @param filepath is a filepath to write to
#'
#' @return \code{IO} nothing
#'
write_to_spreadsheet <- function(lab_book, filepath) {
    for (d in lab_book) {
        xlsx::write.xlsx(x = d$df, file = filepath, sheetName = d$name)
    }
}

