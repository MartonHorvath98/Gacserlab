#! bin/bash

# Take user arguments: input-, config- and output directories
while getopts i:c:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        c) config_file=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

# Read in files from the input directory
if [ -d "$input_dir" ]
then
    # Create the output directory if it does not exist
    if [ ! -d "$output_dir" ]
    then
        mkdir -p $output_dir
    fi

    # Create the samples.tsv file if it does not exist
    if [ ! -f "${config_file}" ]
    then    
        # Loop through the files in the input directory
        #ls ${input_dir}*.fq.gz | sed 's/_1.fq.gz//g' | sed 's/_2.fq.gz//g' |\
        #sed 'n;d' > ${output_dir}${config_file}
        ls ${input_dir}*.fq | sed 's/_raw.fq//g' > ${output_dir}${config_file}
        # remove path from the file names
        sed -i "s|${input_dir}||g" ${output_dir}${config_file}
    fi


else
    echo "Input directory does not exist"
fi