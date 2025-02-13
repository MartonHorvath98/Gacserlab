#!/bin/bash

# This script is used to run the QC pipeline for the RNA-seq data
while getopts i:j:o: flag
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
    echo "Trimming adapters with Trim Galore on sample: ${sample}"
    fw_read=${input_dir}${sample}_1.fq.gz
    rv_read=${input_dir}${sample}_2.fq.gz
    
    if [ ! -f ${out_dir}${sample}_1_val_1.fastq.gz ]
    then
	trim_galore ${fw_read} ${rv_read} -j ${cores} -q 20 -e 0.1 --length 36 --paired \
    --illumina --clip_R1 10 --clip_R2 10 --output_dir ${out_dir}
    fi
done

mkdir-p ${trim_report}
source bin/QC.sh -i ${out_dir} -o ${trim_report}

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit

