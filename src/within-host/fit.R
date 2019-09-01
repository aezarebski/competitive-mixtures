#' Variables provided by the calling environment
#'
#' DUMP_FILE
#' MODEL_FILE
#' IC_FILE
#' OUT_FILE
#'

library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)

source(DUMP_FILE)


#' Return an extended array with the new values equal to the mean.
#'
#' @param a the array
#' @param n the length to extend the first dimension by
extended_by_mean <- function(a, n) {
    dimensionality <- length(dim(a))
    if (dimensionality == 1) {
        a_mean <- mean(a)
        return(as.array(c(a, rep(a_mean, n))))
    } else if (dimensionality == 2) {
        result <- colMeans(a) %>%
            rep(n) %>%
            array(c(ncol(a),n)) %>%
            t() %>%
            {rbind(a, .)}
        return(result)
    } else {
        stop()
    }
}


#' Return an initial condition list based on a stanfit.
#'
#' @param fit is the result of a previous optimisation
#' @param num_pure_ferrets is the number of pure ferrets
#' @param num_mix is the number of mixture ferrets
ic_from_fit <- function(fit, num_pure_ferrets, num_mix) {
    .num_pure_ferrets <- length(fit$par$VRNA_0_pure)
    if (.num_pure_ferrets == num_pure_ferrets) {
        .ic_vrna_0_pure <- fit$par$VRNA_0_pure
        .ic_vtcid_0_pure <- fit$par$VTCID_0_pure
    } else if (.num_pure_ferrets < num_pure_ferrets) {
        .num_missing_pure <- num_pure_ferrets - .num_pure_ferrets
        .ic_vrna_0_pure <- extended_by_mean(fit$par$VRNA_0_pure, .num_missing_pure)
        .ic_vtcid_0_pure <- extended_by_mean(fit$par$VTCID_0_pure, .num_missing_pure)
    } else {
        cat("Not implemented when reducing the number of ferrets\n")
        stop()
    }
    
    .num_mix <- nrow(fit$par$VRNA_0_mix)
    if (.num_mix == num_mix) {
        .ic_vrna_0_mix <- fit$par$VRNA_0_mix
        .ic_vtcid_0_mix <- fit$par$VTCID_0_mix
    } else if (.num_mix < num_mix) {
        .num_missing_mix <- num_mix - .num_mix
        if (.num_mix > 0) {
            .ic_vrna_0_mix <- extended_by_mean(fit$par$VRNA_0_mix, .num_missing_mix)
            .ic_vtcid_0_mix <- extended_by_mean(fit$par$VTCID_0_mix, .num_missing_mix)
        } else {
            .helper <- function(x) {
                #' Return a matrix with mean(x) as values.
                array(rep(mean(x), num_mix * 2), c(num_mix, 2))
            }
            .ic_vrna_0_mix <- .helper(fit$par$VRNA_0_pure)
            .ic_vtcid_0_mix <- .helper(fit$par$VTCID_0_pure)
        }
    } else {
        cat("Not implemented when reducing the number of ferrets\n")
        stop()
    }
    
    return(
        list(
            VRNA_0_pure = .ic_vrna_0_pure,
            VTCID_0_pure = .ic_vtcid_0_pure, 
            VRNA_0_mix = .ic_vrna_0_mix,
            VTCID_0_mix = .ic_vtcid_0_mix, 
            theta = fit$par$theta
        )
    )
}


petrie_stanmodel <- stan_model(MODEL_FILE)
if (IC_FILE != "NA") {

    if (!file.exists(IC_FILE)) {
        cat(sprintf("\n\n\nCannot find initial condition file: %s\n\n\n", IC_FILE))
        stop()
    }

    initial_opt_condition <- IC_FILE %>%
        readRDS() %>%
        ic_from_fit(num_pure_wild + num_pure_mutant, num_mix)
} else {
    initial_opt_condition <- 'random'  # Default value.
}
fit <- optimizing(
    petrie_stanmodel,
    init = initial_opt_condition,
    verbose = TRUE,
    as_vector = FALSE,
    hessian = TRUE
)

if (!dir.exists(dirname(OUT_FILE))) {
    cat(sprintf("\n\n\nCannot find directory: %s\n\n\n", dirname(file)))
    stop()
}


saveRDS(fit, file = OUT_FILE)
