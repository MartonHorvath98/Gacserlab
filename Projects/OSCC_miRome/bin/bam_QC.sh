#!/usr/bin/env bash
# bam_QC.sh
# Marton Horvath
# July, 2024

# This script is used to run QC on the aligned sequence data .bam files
while getopts i:r:g:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        r) ref=${OPTARG};;
        g) config=${OPTARG};;
        o) out_dir=${OPTARG};;
    esac
done

# ------------ Custom download function ------------ #
source ~/gacserlab/bin/make_rRNA.sh ${ref}

# ------------ Set global variables ------------ #
SAMPLES="$(cat ~/gacserlab/config/samples.txt | cut -f1)"
GENES="${ref}Homo_sapiens.GRCh38.bed"
HOUSEKEEPING="${ref}Homo_sapiens.housekeeping.bed"

if [[ ! -s "${GENES}" ]]
then
    agat_convert_sp_gff2bed.pl --gff ${ref}/Homo_sapiens.GRCh38.gff --out ${ref}/Homo_sapiens.GRCh38.bed 
fi

if [[ ! -s "${HOUSEKEEPING}" ]]
then
    # extract transcript IDs
    if mapfile -d\n genes < "${config}/housekeeping.txt";
    then 
        pattern=$(echo ${genes[*]} | tr " " "|")
        cat "${ref}/Homo_sapiens.GRCh38.gtf" | awk ' { if ($3 == "transcript") {print}}' |\
        grep -wE "\b${pattern}\b" | grep -oP "(?<=transcript_id \")\w*" >> "${config}/transcripts.txt"
    fi
    # extract housekeeping genetranscripts from the bed file
    if mapfile -d\n genes < "${config}/transcripts.txt";
    then 
        pattern=$(echo ${genes[*]} | tr " " "|")
        cat ${GENES} | grep -wE "\b${pattern}\b" > "${HOUSEKEEPING}"
    fi
fi
# ------------ Run the workflow ------------ #
## Create miRNA reference
echo "#######################################"
echo " Collecting alignment metrics (Picard) "
echo "#######################################"

mkdir -p ${out_dir} 

for sample in ${SAMPLES}
do
    echo "Collect alignment metrics: ${i}"
    alignment=${input_dir}${sample}.sorted.bam
    
    if [[ ! -f "${input_dir}${sample}/${smple}_metrix.txt" ]]
    then
    picard CollectRnaSeqMetrics -I ${alignment} -O "${input_dir}${sample}/${smple}_metrix.txt" \
    --REF_FLAT "${ref}Homo_sapiens.GRCh38.gtf.refflat" --REFERENCE_SEQUENCE "${ref}Homo_sapiens.GRCh38.fa" \
    --RIBOSOMAL_INTERVALS "${ref}Homo_sapiens.GRCh38.rRNA.interval_list" --STRAND SECOND_READ_TRANSCRIPTION_STRAND --ASSUME_SORTED true
    fi
done

## Create miRNA reference
echo "#######################################"
echo " Calculating coverage and TIN (RSeQC)  "
echo "#######################################"

mkdir -p ${out_dir} 
# calculating the coverage of housekeeping genes
geneBody_coverage.py -r "${HOUSEKEEPING}" -i "${input_dir}"  -o "${output_dir}"

# evaluating RNA integrity at the transcript level
tin.py -r "${HOUSEKEEPING}" -i "${input_dir}"  -o "${output_dir}"

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
