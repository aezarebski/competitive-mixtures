#!/bin/sh

CSV_FILE=in/between-host/Mut-N197data2.csv
DUMP_FILE=out/between-host/mut-d197n.dump
MODEL_FILE=src/between-host/mccaw.stan
FIT_FILE=out/between-host/fit-d197n.rds

echo "Running the inference for the between-host model for D197N"

# If the output directory is missing then exit early with a message.
if [ ! -d "out/between-host" ]; then
    echo "ERROR: Could not find output directory: out/between-host"
    exit 1
fi

echo "Pre-processing data"
Rscript main.R --between_dump --csv $CSV_FILE --dump $DUMP_FILE
echo "Running sampler"
Rscript main.R --between_sample --dump $DUMP_FILE --model $MODEL_FILE --out $FIT_FILE
echo "Post-processing samples"
Rscript main.R --between_view --dump $DUMP_FILE --out $FIT_FILE
echo "Inference finished"
