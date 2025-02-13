#!/bin/bash

# This script is used to run the QC pipeline for the RNA-seq data
while getopts i:o: flag
do
    case "${flag}" in
        i) dir=${OPTARG};;
        o) out_dir=${OPTARG};;
    esac
done

echo "$dir, $out_dir"
# Read in files from the input file or directory using -f or -d flags
if [ -d "$dir" ]
then
    for f in $dir/*.fq.gz; do
        if [ -f "$f" ]; then
            echo "Processing  $(basename $f) file..."
            # Run the QC pipeline
            fastqc $f -t 4 --svg --memory 8192 -o $out_dir
            seqkit stats $f -j 4 -Ta >> ${out_dir}/stats.txt
        fi
    done

    # Clean seqkit output: remove repetitive headers and remove the path from the file name
    stat_file=$out_dir/stats.txt
    awk 'NR==1 || NR%2==0' $stat_file | awk -F'\t' -v OFS='\t' '{sub(".*/", "", $1)} 1' > $out_dir/stats_table.tsv 
    rm $stat_file
else
        echo "The input directory does not exist!"
        exit 1
fi

