#!/bin/bash

# This script is used to run the QC pipeline for the RNA-seq data
while getopts i:j:r:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        j) cores=${OPTARG};;
        r) reference=${OPTARG};;
        o) out_dir=${OPTARG};;
    esac
done

# ------------ Set global variables ------------ #
# Set the number of cores
eval echo "Number of cores: ${cores}"

SAMPLES="$(cat config/samples.txt | cut -f1)"

eval echo "Alignment files are at: ${input_dir}"

# ------------ Run the workflow ------------ #
## Adapter trimming
echo "####################################"
echo "Count features with featureCounts   "
echo "####################################"

## Alignment results directory
mkdir -p ${out_dir}
eval echo "Alignment files: ${out_dir}"

# Adapter trimming
for sample in ${SAMPLES}
do
    echo "Creating alignments for sample: ${sample}"
    alignment=${input_dir}${sample}.sorted.bam
    
    if [ ! -f ${out_dir}${sample}.counts.txt ]
    then
    # Count features with featureCounts
    featureCounts -p --countReadPairs -s 2 -T ${cores} -a ${reference} -o ${out_dir}${sample}.counts.txt ${alignment}
    fi
done

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
