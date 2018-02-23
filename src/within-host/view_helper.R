library(reshape2)
library(dplyr)
library(ggplot2)
library(rstan)


model_solution_as_df <- function(ms) {
    num_ferrets <- dim(ms)[1]
    .f <- function(n) {
        data.frame(
            tcid = log10(ms[n,1,]),
            rna = log10(ms[n,2,]),
            pyro = ms[n,3,],
            ferret_num = n,
            obs_num = 1:length(ms[1,1,])
        )
    }
    lapply(1:num_ferrets, .f) %>%
        {do.call(rbind, .)} %>%
        melt(id.vars = c("ferret_num", "obs_num"))
}


observations_as_df <- function(obs) {
    num_ferrets <- dim(obs)[1]
    num_obs <- dim(obs)[3]
    .f <- function(n) {
        data.frame(
            tcid = obs[n,1,],
            rna = obs[n,2,],
            pyro = obs[n,3,],
            ferret_num = n,
            obs_num = 1:num_obs
        )
    }
    lapply(1:num_ferrets, .f) %>%
        {do.call(rbind, .)} %>%
        melt(id.vars = c("ferret_num", "obs_num"))
}

plot_model_solution <- function(model_solution_array, observations_array) {
    ms_plot_df <- model_solution_as_df(model_solution_array)
    obs_plot_df <- observations_as_df(observations_array)
    ggplot(mapping = aes(x = obs_num, y = value)) +
        geom_line(data = ms_plot_df, colour = "blue") +
        geom_point(data = obs_plot_df, colour = "red") +
        facet_grid(ferret_num~variable) +
        labs(
            x = "Day",
            y = "Observed value"
        )
}


#' Based on the facetted figure create single visualisation for each of the
#' ferrets and save it to a file in the out put directory.
#'
#' @param g_original a ggplot object from plot_model_solution
#' @param out_dir the directory to write the visualisations to.
plot_individual_solutions <- function(g_original, out_dir) {
  gg_data <- ggplot_build(g_original)
  fit_df <- gg_data$data[[1]]
  data_df <- gg_data$data[[2]]
  num_ferrets <- length(unique(data_df$PANEL)) / 3
  
  out_file <- function(ferret_num) {sprintf("%s/ferret_%d.pdf", out_dir, ferret_num)}
  panel_mask <- function(df, ferret_num) {ceiling(as.integer(df$PANEL) / 3) == ferret_num}
  
  for (ferret_num in 1:num_ferrets) {
    mod_frame <- fit_df[panel_mask(fit_df, ferret_num), c("x", "y", "PANEL")]
    obs_frame <- data_df[panel_mask(data_df, ferret_num), c("x", "y", "PANEL")]
    plot_df <- bind_rows(
        mutate(mod_frame, type = "model"),
        mutate(obs_frame, type = "data")) %>%
        rename(variable = PANEL)
    g <- ggplot(mapping = aes(x = x, y = y)) +
        geom_line(data = filter(plot_df, type == "model")) +
        geom_point(data = filter(plot_df, type == "data")) +
        facet_wrap(~variable, scale = "free_y")
    output_file <- out_file(ferret_num)
    ggsave(output_file, g)
  }
}


verbose_ggsave <- function(g, fn) {
    cat(sprintf("Writing to file: %s\n", fn))
    ggsave(fn, g)
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

    

#' Return a vector of ratios of R0 in logarithms.
#'
#' @param constr_samples is a matrix of constrained parameter samples in their
#' natural scale. i.e., not constrained and not logged.
ln_r0_ratio_samples <- function(constr_samples) {
    ln_theta <- log(constr_samples)
    ln_r0_a <- ln_theta[,7] + ln_theta[,1] - (ln_theta[,9] + ln_theta[,5])
    ln_r0_b <- ln_theta[,8] + ln_theta[,2] - (ln_theta[,10] + ln_theta[,6])
    return(ln_r0_a - ln_r0_b)
}


#' Return just the `theta` component of the constrained sample lists as a matrix
#' where each row of the matrix is a sample.
#'
#' @param constr_samples is a list of constrained parameter lists
just_theta <- function(constr_samples) {
    t(sapply(constr_samples, function(cs) {cs[["theta"]]}))
}


#' Return a list of constrained parameter sample lists.
#'
#' @param uncon_samples is a matrix of unconstrained parameter samples
#' @param sf_obj is a `stanfit` object for the petrie model
constrained_samples <- function(uncon_samples, sf_obj) {
    transform_samples <- function(uc_pars) {
        constrain_pars(sf_obj, uc_pars)
    }
    apply(uncon_samples, 2, transform_samples)
}


#' Return a matrix of unconstrained parameter samples.
#'
#' @param fit is an instance of the returned value from `rstan::optimizing`.
#' @param num_samples is the number of samples to generate
#' @param sf_obj is a `stanfit` object for the petrie model
unconstrained_samples <- function(fit, num_samples, sf_obj) {
    uc_mean <- unconstrain_pars(sf_obj, fit$par)
    cat(sprintf("\n\tThe length of the mean vector is %d\n\n", length(uc_mean)))
    uc_hess <- fit$hessian
    cat("\n\tThe dimensions of the Hessian matrix are:\n\n")
    print(dim(uc_hess));
    uc_A <- t(chol(chol2inv(chol(-uc_hess))))
    theta_len <- length(uc_mean)
    print(num_samples);
    z_samples <- rnorm(n = theta_len * num_samples) %>%
        matrix(nrow = theta_len, ncol = num_samples)
    uc_mean + uc_A %*% z_samples
}


#' Return a vector of R_0 ratios from an approximate posterior.
#'
#' @param fit is an instance of the returned value from `rstan::optimizing`.
#' @param params is a list of implementation specifics
r0_ratios <- function(fit, params) {
    fit %>%
        unconstrained_samples(params$num_samples, params$sf_obj) %>%
        constrained_samples(params$sf_obj) %>%
        just_theta() %>%
        ln_r0_ratio_samples() %>%
        exp()
}



#' Return a visualisation of the samples of the ratio of the reproduction numbers.
#'
#' @param r0_ratios is a numeric vector of the ratios in the natural parameters.
plot_posterior_r0_ratios <- function(r0_ratios) {
    est <- quantile(r0_ratios, probs=c(0.025, 0.5, 0.975))
    
    num_geq_1 <- sum(r0_ratios >= 1)
    num_samples <- length(r0_ratios)
    posterior_dist_of_prob_est <- function(p) {
        qbeta(p, shape1 = num_samples - num_geq_1, shape2 = num_geq_1)
    }
    prob_est <- posterior_dist_of_prob_est(c(0.025, 0.5, 0.975))
    prob_est_label <- sprintf("P(ratio < 1) = \n%.02f (%.02f, %.02f)", prob_est[1], prob_est[2], prob_est[3])
    bar <- paste(rep("=", 80), collapse="")
    cat(sprintf("\n%s\n\n", bar))
    cat("\tquantile(r0_ratios, probs=c(0.025, 0.5, 0.975))\n")
    print(est)
    cat("\tposterior_dist_of_prob_est(c(0.025, 0.5, 0.975))\n")
    print(prob_est)
    cat(sprintf("\n%s\n\n", bar))
    
    plot_df <- data.frame(x = r0_ratios)
    ggplot(plot_df, aes(x = x)) +
        geom_histogram() +
        annotate("text", x = 0.75, y = 75, label = prob_est_label) +
        geom_vline(xintercept = 1) +
        labs(x = "R0 ratio") +
        ggtitle(sprintf("R0 ratio and 95%% CI: %.02f ( %.02f, %.02f )", est[2], est[1], est[3]))
}
