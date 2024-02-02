#!/bin/bash

# Grabs a copy of the hg38 reference genome

cd dataPublic
mkdir -p hg38
cd hg38

# Store the split components in a sub-directory
mkdir -p split

url="ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
file_gz="hg38.fna.gz"
file="hg38.fna"

wget -O ${file_gz} ${url}
gunzip ${file_gz}

# Modified from the following:
# https://crashcourse.housegordon.org/split-fasta-files.html

# Split based on '>' character
csplit -s -z ${file} '/>/' '{*}'

# The (temp) output files are labeled as xx00, xx01, etc.
for i in xx*; 
do
    # Get the chromosome ID
    n=$(sed 's/>// ; s/ .*// ; 1q' "$i")
    
    # Rename appropriately
    mv "$i" "split/$n.fa"
done