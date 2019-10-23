#!/usr/bin/env bash

XLSX_FILE=in/within-host/d197n.xlsx
DUMP_FILE_R3=out/within-host/d197n-relax-3.dump
DUMP_FILE_R2=out/within-host/d197n-relax-2.dump
DUMP_FILE_R1=out/within-host/d197n-relax-1.dump
MODEL_FILE=src/within-host/petrie.stan
OUTPUT_FILE_R3=out/within-host/fit-d197n-relax-3.rds
OUTPUT_FILE_R2=out/within-host/fit-d197n-relax-2.rds
OUTPUT_FILE_R1=out/within-host/fit-d197n-relax-1.rds
OUTPUT_DIR=out/within-host

# If the output directory is missing then exit early with a message.
if [ ! -d $OUTPUT_DIR ]; then
    echo "ERROR: Could not find output directory: $OUTPUT_DIR"
    exit 1
fi

echo "Running the demonstration code for the within-host model"
echo "Pre-processing data"
Rscript main.R --within_dump --xlsx $XLSX_FILE --dump $DUMP_FILE_R3 --relax 3
Rscript main.R --within_dump --xlsx $XLSX_FILE --dump $DUMP_FILE_R2 --relax 2
Rscript main.R --within_dump --xlsx $XLSX_FILE --dump $DUMP_FILE_R1 --relax 1
echo "Running the fit"
Rscript main.R --within_fit --dump $DUMP_FILE_R3 --model $MODEL_FILE --out $OUTPUT_FILE_R3
Rscript main.R --within_fit --dump $DUMP_FILE_R2 --model $MODEL_FILE --out $OUTPUT_FILE_R2 --ic $OUTPUT_FILE_R3
Rscript main.R --within_fit --dump $DUMP_FILE_R1 --model $MODEL_FILE --out $OUTPUT_FILE_R1 --ic $OUTPUT_FILE_R2
echo "Post-processing the results"
Rscript main.R --within_view --dump $DUMP_FILE_R1 --model $MODEL_FILE --out $OUTPUT_FILE_R1 --out_dir $OUTPUT_DIR
echo "Demonstation finished"
