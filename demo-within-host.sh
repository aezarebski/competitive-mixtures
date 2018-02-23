#!/usr/bin/env bash

XLSX_FILE=in/demo.xlsx
DUMP_FILE=out/demo/demo.dump
MODEL_FILE=src/within-host/petrie.stan
OUTPUT_FILE=out/demo/fit.rds
OUTPUT_DIR=out/demo
echo "Running the demonstration code for the within-host model"
echo "Pre-processing data"
Rscript main.R --within_dump --xlsx $XLSX_FILE --dump $DUMP_FILE --relax 1
echo "Running the fit"
Rscript main.R --within_fit --dump $DUMP_FILE --model $MODEL_FILE --out $OUTPUT_FILE
echo "Post-processing the results"
Rscript main.R --within_view --dump $DUMP_FILE --model $MODEL_FILE --out $OUTPUT_FILE --out_dir $OUTPUT_DIR
echo "Demonstation finished"
