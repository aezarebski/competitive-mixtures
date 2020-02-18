#' Variables provided by the calling environment
#'
#' DUMP_FILE
#' FIT_FILE
#' OUT_DIR
#' MODEL_FILE
#' POSTERIOR_VRNA_0_SAMPLE_FILE
#' POSTERIOR_MODEL_SOL_SAMPLE_FILE
#'
#'

#' TODO These need to be handled properly!
POSTERIOR_VRNA_0_SAMPLE_FILE <- "demo-vrna.csv"
POSTERIOR_MODEL_SOL_SAMPLE_FILE <- "demo-model.csv"

library(purrr)
library(dplyr)


source("src/within-host/view_helper.R")

#' Visualise the model fit
obs_array <- readRDS(gsub(".dump", "obs_array.rds", DUMP_FILE))
fit <- readRDS(FIT_FILE)
sol <- fit$par$sol
g_model_fit <- plot_model_solution(sol, obs_array)
fit_output_file <- gsub(".rds", ".pdf", gsub("fit", "model-fit", FIT_FILE))
verbose_ggsave(g_model_fit, fit_output_file)

#' Visualise the model fit for each subject
single_fits <- plot_model_solutions_as_list(sol, obs_array)
for (ix in 1:length(single_fits)) {
    verbose_ggsave(single_fits[[ix]]$fig, gsub(".rds", sprintf("-model-ferret-%d.pdf", single_fits[[ix]]$fig_num), FIT_FILE))
}

#' Visualise the posterior distribution on the R0-ratio
source("src/within-host/dumper_helper.R")
source(DUMP_FILE)
params <- list(
    num_samples = 1E2,
    sf_obj = stan(MODEL_FILE, data = get_dump_vars(), chains = 0, iter = 0)
)
r0_samples <- r0_ratios(fit, params)
g_r0_vis <- plot_posterior_r0_ratios(r0_samples)
r0_vis_output_file <- gsub(".rds", ".pdf", gsub("fit", "approx_r0_posterior", FIT_FILE))
verbose_ggsave(g_r0_vis, r0_vis_output_file)

#' Record a copy of some parameter samples from the approximate prior distribution.
post_samples <- posterior_samples(fit, params)

vrna_0_sample_as_df <- function(x, ix) {
    x <- x$VRNA_0_mix %>%
        as.data.frame %>%
        rename(VRNA_0_A = V1, VRNA_0_B = V2)
    x$mix_ferret_num <- 1:nrow(x)
    x$ix <- ix
    return(x)
}
vrna_0_posterior_df <- map2(.x = post_samples, .y = 1:length(post_samples), .f = vrna_0_sample_as_df) %>% bind_rows()
write.table(x = vrna_0_posterior_df, file = POSTERIOR_VRNA_0_SAMPLE_FILE, quote = FALSE, sep = ",", row.names = FALSE)

model_sol_sample_as_df <- function(x, ix) {
    x <- x$model_sol %>%
        reshape2::melt(varnames = c("ferret", "variable", "day")) %>%
        dplyr::mutate(ix = ix)
    return(x)
}
model_sol_posterior_df <- map2(.x = post_samples, .y = 1:length(post_samples), .f = model_sol_sample_as_df) %>% bind_rows()
write.table(x = model_sol_posterior_df, file = POSTERIOR_MODEL_SOL_SAMPLE_FILE, quote = FALSE, sep = ",", row.names = FALSE)



#' Return a list of lists of plots showing the observations and a posterior sample of trajectories.
#'
#' @param obs_df data.frame of melted observation dump array.
#' @param post_df data.frame of posterior samples of the model solution.
#'
#' @examples
#' obs_df <- readRDS("out/demo/demoobs_array.rds") %>% reshape2::melt(varnames = c("ferret", "variable", "day"))
#' post_df <- read.table(file = "demo-model.csv", header = TRUE, sep = ",")
#' all_fits <- all_ferret_fits(obs_df, post_df)
#' print(all_fits[[2]]$pyro)
#'
all_ferret_fits <- function(obs_df, post_df) {

    ferret_plot <- function(f_id) {
        tcid50_and_realtime_ggplot <- ggplot() +
            geom_line(data = dplyr::filter(post_df, variable != 3, ferret == f_id),
                      mapping = aes(x = day, y = log(value), group = ix), colour = "blue", alpha = 0.1) +
            geom_point(data = dplyr::filter(obs_df, variable != 3, ferret == f_id),
                       mapping = aes(x = day, y = value), colour = "red") +
            facet_wrap(~variable) +
            ggtitle(sprintf("TCID50 and RealTime data and posterior fit for ferret %d", f_id))

        # TODO Fix the way that missing data is visualised!

        pyro_ggplot <- ggplot() +
            geom_line(data = dplyr::filter(post_df, variable == 3, ferret == f_id),
                      mapping = aes(x = day, y = value, group = ix), colour = "blue", alpha = 0.1) +
            geom_point(data = dplyr::filter(obs_df, variable == 3, ferret == f_id),
                       mapping = aes(x = day, y = value), colour = "red") +
            ggtitle(sprintf("Pyrosequencing data and posterior fit for ferret %d", f_id))

        list(tcid50_and_realtime = tcid50_and_realtime_ggplot,
             pyro = pyro_ggplot,
             ferret_id = f_id)
    }

    lapply(unique(obs_df$ferret), ferret_plot)
}


obs_df <- readRDS(gsub(".dump", "obs_array.rds", DUMP_FILE)) %>% reshape2::melt(varnames = c("ferret", "variable", "day"))

post_df <- read.table(file = POSTERIOR_MODEL_SOL_SAMPLE_FILE, header = TRUE, sep = ",")

cat("\nSaving all the model fits to all-fits.rds\n")
saveRDS(object = all_ferret_fits(obs_df, post_df), file = "all-fits.rds")
