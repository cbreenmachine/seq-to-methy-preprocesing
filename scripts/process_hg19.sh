#!/bin/bash

# Follow these instruction
# https://www.gungorbudak.com/blog/2014/04/13/download-human-reference-genome-hg19/
url="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz"

cd dataPublic
mkdir -p hg19
cd hg19

wget -O hg19.raw.fna.gz ${url}
gunzip hg19.raw.fna.gz

# cat hg19.raw.fna \
#     | sed -n '/>NC./,/>./{/>./d;p}' \
#     | sed 's/.*chromosome \(X\|Y\|[0-9]*\).*/>chr\1/' > hg19.fna