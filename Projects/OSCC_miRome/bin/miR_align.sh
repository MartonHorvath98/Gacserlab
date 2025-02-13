#!/bin/bash

# This script is used to run the QC pipeline for the RNA-seq data
while getopts i:r:j:b:m:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        r) genome=${OPTARG};;
        j) cores=${OPTARG};;
        b) bowtie=${OPTARG};;
        m) mirbase=${OPTARG};;
        o) out_dir=${OPTARG};;
    esac
done

# ------------ Custom download function ------------ #
get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o $2 $1
		return $?
	else
		wget $1 -P $2
		return $?
	fi
}
# ------------ Set global variables ------------ #
SAMPLES="$(cat ~/gacserlab/config/samples.txt | cut -f1)"
MIRBASE="https://www.mirbase.org/download/CURRENT/"
MATURE="mature.fa"
HAIRPIN="hairpin.fa"

# Set the number of cores
eval echo "Number of cores: ${cores}"
eval echo "Trimmmed read are at: ${input_dir}"

# ------------ Run the workflow ------------ #
## Create miRNA reference
if [[ ! -s "${mirbase}mature_ref.fa" ]]
then
    # Download and parse the miRBase mature reference
    get ${MIRBASE}${MATURE} ${mirbase} || (echo "Error getting ${MATURE}" && exit 1)
    # Extract the human mature miRNAs
    cat ${mirbase}${MATURE} | sed -e 's/&gt//g' | tr ";" "\n" | grep -E "\bhsa\b" |\
    awk -F '<br>' '{split($1, a, " "); printf(">%s\n%s\n", a[1], $2)}' > "${mirbase}${MATURE%.fa}_ref.fa"
    # Extract other mature miRNAs (mouse, and chimpanzee)
    cat ${mirbase}${MATURE} | sed -e 's/&gt//g' | tr ";" "\n" | grep -E "\bmmu|ptr\b" |\
    awk -F '<br>' '{split($1, a, " "); printf(">%s\n%s\n", a[1], $2)}' > "${mirbase}${MATURE%.fa}_other.fa"
    rm ${mirbase}${MATURE}

    # Download and parse the miRBase hairpin reference
    get ${MIRBASE}${HAIRPIN} ${mirbase} || (echo "Error getting ${HAIRPIN}" && exit 1)
    # Extract the human hairpin miRNAs
    cat ${mirbase}${HAIRPIN} | sed -e 's/&gt//g' | tr ";" "\n" | grep -E "\bhsa\b" |\
    awk -F '<br>' '{split($1, a, " "); printf(">%s\n%s\n", a[1], $2)}' > "${mirbase}${HAIRPIN%.fa}_ref.fa"
    rm ${mirbase}${HAIRPIN}
fi

## Create bowtieindex for the reference genome
if [[ ! -s "${bowtie}.1.ebwt" ]]
then
    echo "Creating bowtie index for the reference genome"
    awk '/^>/ {if ($1 !~ /^>[KG]/) sub(/^>/, "&chr"); print $1; next} {print}' ${genome} > ${mirbase}/hg38.fa 
    bowtie-build --seed 100 --threads 16 ${mirbase}/hg38.fa ${bowtie}
fi

echo "####################################"
echo "Align reads with miRdeep2           "
echo "####################################"

mkdir -p ${out_dir}

for sample in ${SAMPLES}
do
    echo "Creating alignments for sample: ${sample}"
    read=${input_dir}${sample}_trimmed.fq
    
    # if [ ! -f ${out_dir}${sample}/${sample}.arf ]
    # then
    mkdir -p ${out_dir}${sample} && cd ${out_dir}${sample}
    ## Collapse identical miRNA reads
    echo "Collapsing identical miRNA reads"
    mapper.pl ${read} -e -h -i -j -m -u -n -o ${cores} -p ${bowtie} -s ${out_dir}${sample}/${sample}_collapsed.fa \
    -t ${out_dir}${sample}/${sample}.arf 2> ${out_dir}${sample}/${sample}_mapper.log

    # Quantify the miRNA expression
    # echo "Quantifying miRNA expression"
    # miRDeep2.pl ${out_dir}${sample}/${sample}_collapsed.fa ${mirbase}/hg38.fa ${out_dir}${sample}/${sample}.arf \
    # ${mirbase}${MATURE%.fa}_ref.fa ${mirbase}${MATURE%.fa}_other.fa ${mirbase}${HAIRPIN%.fa}_ref.fa \
    # -t Human -v -c -P 2> ${out_dir}${sample}/${sample}_align.log

    cd ${out_dir}
    # fi
done

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
