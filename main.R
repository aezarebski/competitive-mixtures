library(optparse)
set.seed(1)


#' Parse the command line arguments to determine what to do.
get_args <- function() {
    opts <- list(
        make_option(
            c("--within_dump"),
            action = "store_true",
            type = "logical",
            default = FALSE,
            help = "Dump the data for the within-host model."
        ),
        make_option(
            c("--xlsx"),
            type = "character",
            default = "NA",
            help = "XLSX file containing the raw measurements."
        ),
        make_option(
            c("--dump"),
            type = "character",
            default = "NA",
            help = "DUMP file containing the preprocessed data."
        ),
        make_option(
            c("--relax"),
            type = "numeric",
            default = 1,
            help = "Relaxation factor (defaults to 1 for no relaxation)."
        ),
        make_option(
            c("--within_fit"),
            action = "store_true",
            type = "logical",
            default = FALSE,
            help = "Fit the within-host model."
        ),
        make_option(
            c("--model"),
            type = "character",
            default = "NA",
            help = "The STAN model file to use."
        ),
        make_option(
            c("--out"),
            type = "character",
            default = "NA",
            help = "The RDS file to write the fit to."
        ),
        make_option(
            c("--within_view"),
            action = "store_true",
            default = FALSE,
            help = "Post-process the fit from the within-host model."
        ),
        make_option(
            c("--out_dir"),
            type = "character",
            default = "NA",
            help = "The directory to write output figures to."
        )
    )
    return(parse_args(OptionParser(option_list = opts)))
}


main <- function() {
    args <- get_args()
    if (args$within_dump) {
        dump_env <- new.env()
        dump_env$DATA_FILE <- args$xlsx
        dump_env$DUMP_FILE <- args$dump
        dump_env$RELAX_FACT <- args$relax
        source("src/within-host/dumper.R", local = dump_env)
    } else if (args$within_fit) {
        fit_env <- new.env()
        fit_env$DUMP_FILE <- args$dump
        fit_env$MODEL_FILE <- args$model
        fit_env$OUT_FILE <- args$out
        fit_env$IC_FILE <- "NA"
        source("src/within-host/fit.R", local = fit_env)
    } else if (args$within_view) {
        view_env <- new.env()
        view_env$DUMP_FILE <- args$dump
        view_env$MODEL_FILE <- args$model
        view_env$FIT_FILE <- args$out
        view_env$OUT_DIR <- args$out_dir
        source("src/within-host/view.R", local = view_env)
    }
}


main()

