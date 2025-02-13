#!/bin/bash

# Downloads sequences from Candida Genome Database

# Take user arguments: input directory
while getopts d:j:b: flag
do
    case "${flag}" in
        d) dir=${OPTARG};;
		j) threads=${OPTARG};;
    esac
done

# Set up the URLs
ENSEMBL_RELEASE=58
FASTA_REF=ftp://ftp.ensemblgenomes.org/pub/fungi/release-${ENSEMBL_RELEASE}/fasta/

# Custom download function
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

# Check command executable paths
## 1. BOWTIE2
BOWTIE2_BUILD_EXE=./bowtie2-build
if [ ! -x "$BOWTIE2_BUILD_EXE" ] ; then
	if ! which bowtie2-build ; then
		echo "Could not find bowtie2-build in current directory or in PATH"
		exit 1
	else
		BOWTIE2_BUILD_EXE=`which bowtie2-build`
	fi
fi

# Download the Candida albicans reference genome fasta file (unless it already exists)
ALBICANS_GENOME=Candida_albicans.GCA000182965v3.dna.toplevel.fa.gz
if [ ! -f "${dir}/${ALBICANS_GENOME%.GCA000182965v3.dna.toplevel.fa.gz}.fa" ] ; then
	get ${FASTA_REF}/candida_albicans/dna/$ALBICANS_GENOME ${dir} || (echo "Error getting $ALBICANS_GENOME" && exit 1)
	gunzip ${dir}/$ALBICANS_GENOME || (echo "Error unzipping $ALBICANS_GENOME" && exit 1)
	mv ${dir}/${ALBICANS_GENOME%.gz} ${dir}/${ALBICANS_GENOME%.GCA000182965v3.dna.toplevel.fa.gz}.fa
fi
ALBICANS_GENOME=${dir}/${ALBICANS_GENOME%.GCA000182965v3.dna.toplevel.fa.g}.fa

# Download the Candida parapsilosis reference genome fasta file (unless it already exists)
PARAPSILOSIS_GENOME=Candida_parapsilosis.GCA000182765v2.dna.toplevel.fa.gz
if [ ! -f "${dir}/${PARAPSILOSIS_GENOME%.GCA000182765v2.dna.toplevel.fa.gz}.cds.fa" ] ; then
	get ${ENSEMBL_TRANSCRIPTOME}/candida_parapsilosis/dna/$PARAPSILOSIS_GENOME ${dir} || (echo "Error getting $PARAPSILOSIS_GENOME" && exit 1)
	gunzip ${dir}/$PARAPSILOSIS_GENOME || (echo "Error unzipping $PARAPSILOSIS_GENOME" && exit 1)
	mv ${dir}/${PARAPSILOSIS_GENOME%.gz} ${dir}/${PARAPSILOSIS_GENOME%.GCA000182765v2.dna.toplevel.fa.gz}.fa
fi
PARAPSILOSIS_GENOME=${dir}/${PARAPSILOSIS_GENOME%.GCA000182765v2.dna.toplevel.fa.gz}.fa

echo "##########################################################"
echo " Building BOWTIE2 index                                  #"
echo "##########################################################"

BOWTIE2="${BOWTIE2_BUILD_EXE} --threads ${threads} --seed 100 -f"
echo Running $BOWTIE2 $ALBICANS_GENOME "C_albicans"
if $BOWTIE2 $ALBICANS_GENOME "${dir}/C_albicans" ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi

if $BOWTIE2 $PARAPSILOSIS_GENOME "${dir}/C_parapsilosis" ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi


echo "##########################################################"
echo "# Indexing complete                                      #"
echo "##########################################################"
