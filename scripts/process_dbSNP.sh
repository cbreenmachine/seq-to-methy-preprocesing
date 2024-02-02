#!/bin/bash
# Download and format dbSNP data.
# Specifically, (i) rename chromosomes from 1 to chr1, (ii) index where needed


cd dataPublic/
mkdir -p dbSNP
cd dbSNP

url="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz"
ofile="common_all_20180418.vcf.gz"

map="N_to_chrN.txt"
ofile_renamed="common_all_20180418.renamed.vcf.gz"

for ii in {1..22}; do 
    echo "${ii} chr${ii}" >> ${map}
done

# Download without the sub-directories
wget -O ${ofile} ${url}
bcftools index -f ${ofile}

bcftools annotate --rename-chrs ${map} ${ofile} | bgzip > ${ofile_renamed}
bcftools index -f ${ofile_renamed}
