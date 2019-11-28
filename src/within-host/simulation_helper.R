#' Helper functions for simulating ferret competitive mixture experiments.
library(testthat)


#' \strong{IO}: Observe the ferret
#'
#' @param hidden_state is the internal state of the ferret from \code{ferret}
#' @param observer is the parameters of the observation model
#'
#' @return \code{IO} observations
#'
observe_ferret <- function(hidden_state, observer) {
    tcid_obs <- rnorm(n = nrow(hidden_state),
                      mean = log10(hidden_state[,"TCID_comb"]),
                      sd = observer$sigma_tcid)
    rna_obs <- rnorm(n = nrow(hidden_state),
                     mean = log10(hidden_state[,"RNA_comb"]),
                     sd = observer$sigma_rna)
    pyro_obs <- rnorm(n = nrow(hidden_state),
                      mean = hidden_state[,"ratio_B"],
                      sd = observer$sigma_pyro)
    f_obs <- matrix(c(tcid_obs,rna_obs,pyro_obs), ncol=3)
    colnames(f_obs) <- c("tcid", "rna", "pyro")
    return(f_obs)
}

.tiv_initial_state <- function(params) {
    c(T = params$num_target_cells,
      LA = 0,
      IA = 0,
      VTCIDA = params$v_tcid_a_0,
      VRNAA = params$rho_0 * params$v_tcid_a_0,
      LB = 0,
      IB = 0,
      VTCIDB = params$v_tcid_b_0,
      VRNAB = params$rho_0 * params$v_tcid_b_0)
}

.tiv_differential <- function(params) {
    function(t, y, parms, ...) {
        list(c(-parms$betaA * y["T"] * y["VTCIDA"] - parms$betaB * y["T"] * y["VTCIDB"],
               parms$betaA * y["T"] * y["VTCIDA"] - parms$kA * y["LA"],
               parms$kA * y["LA"] - parms$deltaA * y["IA"],
               parms$pA * y["IA"] - parms$chA * y["VTCIDA"] - parms$dINFA * y["VTCIDA"],
               parms$xiA * parms$pA * y["IA"] - parms$chA * y["VRNAA"],
               parms$betaB * y["T"] * y["VTCIDB"] - parms$kB * y["LB"],
               parms$kB * y["LB"] - parms$deltaB * y["IB"],
               parms$pB * y["IB"] - parms$chB * y["VTCIDB"] - parms$dINFB * y["VTCIDB"],
               parms$xiB * parms$pB * y["IB"] - parms$chB * y["VRNAB"]))
    }
}

.observation_times <- function(protocol) {
    NULL
}

.tiv_params <- function(params) {
    list(
        betaA = params$beta_a,
         kA = params$k_a,
         deltaA = params$delta_a,
         pA = params$p_a,
         chA = params$c_h_a,
         dINFA = params$d_inf_a,
         xiA = params$xi_a,
        betaB = params$beta_b,
        kB = params$k_b,
        deltaB = params$delta_b,
        pB = params$p_b,
        chB = params$c_h_b,
        dINFB = params$d_inf_b,
        xiB = params$xi_b
         )
}

.observation_times <- function(protocol) {
    seq(from = 0, to = protocol$num_observations - 1, by = 1.0)
}


#' The hidden state of the ferret throughout an experiment
#'
#' @param params is the parameters of the ODE model
#' @param protocol is the parameters of the experimental setup
#'
#' @return hidden state of ferret which is a matrix with a row for each observation and a column for each statistic.
#'
ferret <- function(params, protocol) {
    .tiv0 <- .tiv_initial_state(params)
    .diff <- .tiv_differential(params)
    .parms <- .tiv_params(params)
    .obs_times <- .observation_times(protocol)
    expect_that(.diff(.obs_times[1], .tiv0, .parms), is_a("list"))
    expect_that(.diff(.obs_times[1], .tiv0, .parms)[[1]], is_a("numeric"))
    tiv_sol <- deSolve::ode(.tiv0, .obs_times, .diff, .parms)
    expect_that(tiv_sol, is_a("matrix"))
    tcid_comb <- tiv_sol[,"VTCIDA"] + tiv_sol[,"VTCIDB"]
    rna_comb <- tiv_sol[,"VRNAA"] + tiv_sol[,"VRNAB"]
    ratio_b <- tiv_sol[,"VRNAB"] / rna_comb
    fm <- matrix(c(tcid_comb,rna_comb,ratio_b), ncol = 3)
    colnames(fm) <- c("TCID_comb", "RNA_comb", "ratio_B")
    return(fm)
}


#' \strong{IO}: Run the experiment and return the observations.
#
#' @param params is the parameters of the ODE model
#' @param protocol is the parameters of the experimental setup
#' @param observer is the parameters of the observation model
#'
#' @return \code{IO} observations
run_experiment <- function(params, protocol, observer) {
    expect_that(params, has_names(
                            c("num_target_cells",
                              "v_tcid_a_0",
                              "v_tcid_b_0",
                              "rho_0",
                              "beta_a",
                              "k_a",
                              "delta_a",
                              "p_a",
                              "c_h_a",
                              "d_inf_a",
                              "xi_a",
                              "beta_b",
                              "k_b",
                              "delta_b",
                              "p_b",
                              "c_h_b",
                              "d_inf_b",
                              "xi_b")
                          , ignore.order = TRUE))
    expect_that(protocol, has_names(c("num_observations"), ignore.order = TRUE))
    expect_that(observer, has_names(c("sigma_tcid", "sigma_rna", "sigma_pyro"), ignore.order = TRUE))
    hidden_state <- ferret(params, protocol)
    observe_ferret(hidden_state, observer)
}

#' Results of an experiment in a convenient format
#'
#' @param ferrets_and_observations lists of pairs of ferret identifiers and observations
#'
#' @return laboratory book which is a list of lists with a sheet name and the values in that sheet.
#'
laboratory_book <- function(ferrets_and_observations) {
    .f <- function(f_and_o) {
        fid <- f_and_o[[1]]
        obs <- f_and_o[[2]]
        tcid_df <- eval(parse(text = sprintf("data.frame(TCID_%s = obs[,\"tcid\"])", fid)))
        rna_df <- eval(parse(text = sprintf("data.frame(RNA_%s = obs[,\"rna\"])", fid)))
        pyro_df <- eval(parse(text = sprintf("data.frame(PYRO_%s = obs[,\"pyro\"])", fid)))
        list(tcid_df, rna_df, pyro_df)
    }
    .transposed_dfs <- lapply(ferrets_and_observations, .f)

    tcid_df <- do.call(cbind,lapply(.transposed_dfs, function(tdf) tdf[[1]]))
    rna_df <- do.call(cbind,lapply(.transposed_dfs, function(tdf) tdf[[2]]))
    pyro_df <- do.call(cbind,lapply(.transposed_dfs, function(tdf) tdf[[3]]))

    list(
        list(sheet = "TCID50",
              values = do.call(cbind,lapply(.transposed_dfs, function(tdf) tdf[[1]]))),
        list(sheet = "RealTime",
             values = do.call(cbind,lapply(.transposed_dfs, function(tdf) tdf[[2]]))),
        list(sheet = "PYRO",
             values = do.call(cbind,lapply(.transposed_dfs, function(tdf) tdf[[3]])))
         )
}

#' \strong{IO}: Write the results to a spreadsheet
#'
#' @param lab_book is a labbook from \code{laboratory_book}
#' @param filepath is a filepath to write to
#'
#' @return \code{IO} nothing
#'
write_to_spreadsheet <- function(lab_book, filepath) {
    for (ix in 1:3) {
        xlsx::write.xlsx(x = lab_book[[ix]]$values, file = filepath, sheetName = lab_book[[ix]]$sheet, append = TRUE)
    }
}

