#!/bin/bash

# This script is used to run the QC pipeline for the RNA-seq data
while getopts i:j:o:t: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        j) cores=${OPTARG};;
        o) out_dir=${OPTARG};;
        t) trim_report=${OPTARG};;
    esac
done

# ------------ Set global variables ------------ #
# Set the number of cores
eval echo "Number of cores: ${cores}"

SAMPLES="$(cat config/samples.txt | cut -f1)"

eval echo "Raw data at: ${input_dir}"

# ------------ Run the workflow ------------ #
## Adapter trimming
echo "####################################"
echo "Adapter trimming with TrimGalore    "
echo "####################################"

## Create the adapter trimming results directory
mkdir -p ${out_dir}
eval echo "Trimmed reads: ${out_dir}"

# Adapter trimming
for sample in ${SAMPLES}
do
    fw_read=${input_dir}${sample}_raw.fq
    
    if [ ! -f ${out_dir}${sample}_trimmed.fq ]
    then
        echo "Trimming adapters with Trim Galore on sample: ${sample}"

        trim_galore ${fw_read} -j ${cores} --phred33 -q 20 -e 0.1 --length 18 --basename ${sample} \
        --illumina --fastqc_args "--outdir ${trim_report}" --output_dir ${out_dir}
    else 
        echo "Trimmed reads already exist for sample: ${sample}"
    fi
done


# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit