#!/usr/bin/env bash
# make_rRNA.sh
# Kamil Slowikowski
# December 12, 2014
#
# 1. Download chromosome sizes from UCSC if needed.
# 2. Make an interval_list file suitable for CollectRnaSeqMetrics.jar.
#
# Set up the Ensembl release version and base URLs
# ENSEMBL_RELEASE=96
# ENSEMBL_GENES=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/
#
# Picard Tools CollectRnaSeqMetrics.jar:
#
#   https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics

# 1. Chromosome sizes from the UCSC genome browser categories: notes
# ---------------------------
reference=$1
chrom_sizes=Homo_sapiens.GRCh38.chrom_sizes

if [[ ! -s "${reference}${chrom_sizes}" ]]
then
    cut -f1,2 ${reference}Homo_sapiens.GRCh38.fa.fai > ${reference}${chrom_sizes}
fi

# 2. rRNA interval_list file categories: notes
# ------------------------------------------------

# Output file suitable for Picard CollectRnaSeqMetrics.jar.
genes=Homo_sapiens.GRCh38.gtf
rRNA=Homo_sapiens.GRCh38.rRNA.interval_list
 
# Sequence names and lengths. (Must be tab-delimited.)
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:GRCh38"' "${reference}${chrom_sizes}" | \
    grep -v _ \
>> "${reference}${rRNA}"

# Intervals for rRNA transcripts.
grep 'gene_biotype "rRNA"' "${reference}${genes}" | \
    awk '$3 == "transcript"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '
        /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n \
>> "${reference}${rRNA}"   