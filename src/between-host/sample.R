library(rstan)
rstan_options(auto_write = TRUE)


data_list <- readRDS(DUMP_FILE)
fit <- stan(
    MODEL_FILE,
    data = data_list
)
saveRDS(fit, OUT_FILE)
    
