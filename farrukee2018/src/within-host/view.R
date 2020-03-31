#' Variables provided by the calling environment
#'
#' DUMP_FILE
#' FIT_FILE
#' OUT_DIR
#' MODEL_FILE
#'

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

    
