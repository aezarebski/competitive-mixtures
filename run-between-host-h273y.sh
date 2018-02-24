#!/usr/bin/env bash

CSV_FILE=in/between-host/Mut-H273data2.csv
DUMP_FILE=out/between-host/mut-h273y.dump
MODEL_FILE=src/between-host/mccaw.stan
FIT_FILE=out/between-host/fit-h273y.rds

echo "Running the inference for the between-host model for H273Y"
echo "Pre-processing data"
Rscript main.R --between_dump --csv $CSV_FILE --dump $DUMP_FILE
echo "Running sampler"
Rscript main.R --between_sample --dump $DUMP_FILE --model $MODEL_FILE --out $FIT_FILE
echo "Post-processing samples"
Rscript main.R --between_view --dump $DUMP_FILE --out $FIT_FILE
echo "Inference finished"
