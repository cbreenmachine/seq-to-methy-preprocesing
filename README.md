# seq-to-methy-preprocesing
Running the scripts to encode variants, filter by coverage, and standardize covariates.


# Required Files

To run the pipeline, a few starting files are neccessary:
1. Variant call files in `dataRaw/variant-calls/xxx.pass.vcf` are VCFs that passed quality filters
2. Methylation count files in `dataRaw/methylation/xxx.bed` have columns chrom, start, end, methylated, coverage.
3. A chromosome renaming mapping in `dataRaw/dbSNP/N_to_chrN.txt` which has one line per autosome like `1 chr1`. This is used to conver dbSNP data (which used NCBI convention of 1, 2, 3 to label the first three chromosomes) to the system in these VCFs. 


A reference genome and dbSNP data is downloaded automatically.

# Setup

Because snakemake needs to know the output files beforehand, and in this workflow the output files depend on randomization, we need to run the following 

```
python scripts/randomize_samples.py \
        --n_train 185 \
        --n_valid 11 \
        --n_test 185 \
        --ofile "dataDerived/randomization.csv"
```

# Organization

`dataRaw/` contains three sub-directories: `methylation/`, which contains files labeled `{sample_id}.bed`; `variant-calls/`, which contains VCF files similarly named `{sample_id}.pass.vcf`; and `phenotypes/` which contains a CSV with sample covariates such as age, sex, BMI, and disease status. `sample_id`s are three-digit numbers.

`dataDerived/` contains the processed (encoded) data. Files are named as `{sample_id}.{chromosome}.{encoding_scheme}.pt`. Briefly, these files contain tensors 

`scripts/` contains R and Python scripts to munge data. It would be best practice to make this all python, but I'm sticking with what I know for now.


`snakemake  --rerun-triggers mtime -n`