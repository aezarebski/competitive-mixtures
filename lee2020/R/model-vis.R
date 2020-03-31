library(dplyr)
library(reshape2)
library(ggplot2)

save_figures <- TRUE
output_dir <- "doc/img/"

model_output_file <- "foobar.ssv"

if (!file.exists(model_output_file)) {
    stop(sprintf("Could not find model solution file: %s", model_output_file))
}



x <- read.table(model_output_file, header = F, sep = " ")

names(x) <- c("hT",
              "hLA",
              "hLB",
              "hIA",
              "hIB",
              "vTCIDA",
              "vTCIDB",
              "vRNAA",
              "vRNAB")

time_mesh <- seq(from = 0, to = 10, length = nrow(x))




host_plot_df <- x %>% select(starts_with("h")) %>% mutate(time = time_mesh) %>% melt(id.vars = "time")

host_fig <- ggplot(host_plot_df, aes(x = time, y = log(value), colour = variable)) + geom_line()

if (save_figures) {
    ggsave(sprintf("%s/host-fig.pdf", output_dir), plot = host_fig, width = 14.8, height = 10.5, units = "cm")
}



virus_plot_df <- x %>% select(starts_with("v")) %>% mutate(time = time_mesh) %>% melt(id.vars = "time")

virus_fig <- ggplot(virus_plot_df, aes(x = time, y = log(value), colour = variable)) + geom_line()

if (save_figures) {
    ggsave(sprintf("%s/virus-fig.pdf", output_dir), plot = virus_fig, width = 14.8, height = 10.5, units = "cm")
}

