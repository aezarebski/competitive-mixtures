library(optparse)
set.seed(1)


#' Return a logical option for parsing command line arguments.
#'
#' @param flag is a long flag
#' @param help_msg is the documentation for this argument
logical_option <- function(flag, help_msg, ...) {
    make_option(
        flag,
        action = "store_true",
        type = "logical",
        default = FALSE,
        help = help_msg,
        ...
    )
}


#' Return a file option for parsing command line arguments.
#'
#' @param flag is a long flag
#' @param help_msg is the documentation for this argument
file_option <- function(flag, help_msg, ...) {
    make_option(
        flag,
        type = "character",
        default = "NA",
        help = help_msg,
        ...
    )
}


#' Parse the command line arguments to determine what to do.
get_args <- function() {
    opts <- list(
        logical_option(
            "--within_dump",
            "Dump the data for the within-host model."
        ),
        logical_option(
            "--within_fit",
            "Fit the within-host model."
        ),
        logical_option(
            "--within_view",
            "Post-process the fit from the within-host model."
        ),
        logical_option(
            "--between_dump",
            "Dump the data for the between-host model."
        ),
        logical_option(
            "--between_sample",
            "Sample from the between-host model."
        ),
        logical_option(
            "--between_view",
            "Post-process the samples from the between-host model."
        ),
        file_option(
            "--xlsx",
            "XLSX file containing the raw measurements."
        ),
        file_option(
            "--csv",
            "CSV file containing the raw measurements."
        ),
        file_option(
            "--dump",
            "DUMP file containing the preprocessed data."
        ),
        file_option(
            "--model",
            "The STAN model file to use."
        ),
        file_option(
            "--out",
            "The RDS file to write the fit to."
        ),
        file_option(
            "--ic",
            "The RDS file to use as an initial condition."
        ),
        make_option(
            c("--relax"),
            type = "numeric",
            default = 1,
            help = "Relaxation factor (defaults to 1 for no relaxation)."
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
        fit_env$IC_FILE <- args$ic
        source("src/within-host/fit.R", local = fit_env)
    } else if (args$within_view) {
        view_env <- new.env()
        view_env$DUMP_FILE <- args$dump
        view_env$MODEL_FILE <- args$model
        view_env$FIT_FILE <- args$out
        view_env$OUT_DIR <- args$out_dir
        source("src/within-host/view.R", local = view_env)
    } else if (args$between_dump) {
        dump_env <- new.env()
        dump_env$CSV_FILE <- args$csv
        dump_env$DUMP_FILE <- args$dump
        source("src/between-host/dumper.R", local = dump_env)
    } else if (args$between_sample) {
        sample_env <- new.env()
        sample_env$DUMP_FILE <- args$dump
        sample_env$MODEL_FILE <- args$model
        sample_env$OUT_FILE <- args$out
        source("src/between-host/sample.R", local = sample_env)
    } else if (args$between_view) {
        view_env <- new.env()
        view_env$DUMP_FILE <- args$dump
        view_env$OUT_FILE <- args$out
        source("src/between-host/view.R", local = view_env)
    }
}


main()

