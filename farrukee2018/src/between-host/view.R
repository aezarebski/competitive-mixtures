library(rstan)
library(ggplot2)


#' Generate and save diagnostic plots for the stan samples
#'
#' @param fit is a sample object from stan
#' @param output_file_pattern is a `sprintf` type pattern for the names of
#' output files.
make_diagnostic_plots <- function(fit, output_file_pattern) {
    ggsave(
        sprintf(output_file_pattern, "rhat"),
        stan_rhat(fit)
    )
    ggsave(
        sprintf(output_file_pattern, "ess"),
        stan_ess(fit)
    )
    ggsave(
        sprintf(output_file_pattern, "mcse"),
        stan_mcse(fit)
    )
}


#' Return a data frame representing the model fit.
#' 
#' @param fit is a sample object from stan
sample_data_as_plot_df <- function(fit) {
    theta_samples <- as.data.frame(fit)$theta
    exp_neg_theta <- exp(-theta_samples)
    f <- function(theta, d) {
        # Theta is the parameter of the model and d the donor proportion.
        d / (d + (1 - d) * exp(theta))
    }
    summary_stats <- quantile(
        x = theta_samples,
        probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
    )
    d_vals <- seq(1E-6, 1 - 1E-6, length.out = 100)
    data.frame(
        x = d_vals,
        y_min = f(summary_stats[1], d_vals),
        y_low = f(summary_stats[2], d_vals),
        y_med = f(summary_stats[3], d_vals),
        y_hig = f(summary_stats[4], d_vals),
        y_max = f(summary_stats[5], d_vals)
    )
}


#' Return a data frame representing the data the model was fit to.
#'
#' @param data is the data the model was fit to
input_data_as_plot_df <- function(data) {
    data.frame(x = data$d_hat, y = data$r_hat)
}


#' Return a visualisation of the model fit to the data.
#'
#' @param fit is a sample object from stan
#' @param data is the data the model was fit to
model_fit_plot <- function(fit, data) {
    fit_plot_df <- sample_data_as_plot_df(fit)
    data_plot_df <- input_data_as_plot_df(data)
    ggplot() +
        geom_segment(
            mapping = aes(x = 0, y = 0, xend = 1, yend = 1)
        ) +
        geom_ribbon(
            data = fit_plot_df,
            mapping = aes(x = x, ymin = y_min, ymax = y_max),
            alpha = 0.1,
            fill = "blue"
        ) +
        geom_ribbon(
            data = fit_plot_df,
            mapping = aes(x = x, ymin = y_low, ymax = y_hig),
            alpha = 0.1,
            fill = "blue"
        ) +
        geom_line(
            data = fit_plot_df,
            mapping = aes(x = x, y = y_med),
            colour = "blue"
        ) +
        geom_point(
            data = data_plot_df,
            mapping =  aes(x = x, y = y),
            colour = "red"
        ) +
        labs(
            x = "Observed donor proportion",
            y = "Observed recipeint proportion"
        )
}


data_list <- readRDS(DUMP_FILE)
fit <- readRDS(OUT_FILE)

diagnostic_file_pattern <- gsub(".rds", "-diagnostic-%s.pdf", OUT_FILE)
make_diagnostic_plots(fit, diagnostic_file_pattern)

g <- model_fit_plot(fit, data_list)
model_fit_file <- gsub(".rds", "-model-figure.pdf", OUT_FILE)
ggsave(model_fit_file, g)
