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

eval echo "Trimmmed read are at: ${input_dir}"

# ------------ Run the workflow ------------ #
## Adapter trimming
echo "####################################"
echo "Read alginment with Hisat2          "
echo "####################################"

## Alignment results directory
mkdir -p ${out_dir}
eval echo "Alignment files: ${out_dir}"

# Adapter trimming
for sample in ${SAMPLES}
do
    echo "Creating alignments for sample: ${sample}"
    fw_read=${input_dir}${sample}_1_val_1.fq.gz
    rv_read=${input_dir}${sample}_2_val_2.fq.gz
    
    if [ ! -f ${out_dir}${sample}.bam ]
    then
    mkdir -p ${out_dir}${sample}
    # Align reads with Hisat2
    hisat2 --phred33 --dta --non-deterministic --rna-strandness "RF" -p ${cores} \
    --novel-splicesite-outfile "${out_dir}${sample}/${sample}_novel_splicesites.txt" \
    --summary-file "${out_dir}${sample}/stats.txt" --new-summary \
    -x "${reference}" -1 ${fw_read} -2 ${rv_read} |\
    samtools view -h -bS > ${out_dir}${sample}.bam

    # Sort the BAM file
    samtools sort -@ ${cores} -o ${out_dir}${sample}.sorted.bam ${out_dir}${sample}.bam

    # Index the sorted BAM file
    samtools index ${out_dir}${sample}.sorted.bam

    # Remove the unsorted BAM file
    rm ${out_dir}${sample}.bam
    fi
done

mkdir -p ${trim_report}
source bin/QC.sh -i ${out_dir} -o ${trim_report}

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
