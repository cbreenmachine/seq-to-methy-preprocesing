# seq-to-methy-preprocesing
Running the scripts to encode variants, filter by coverage, and standardize covariates.

# Organization

`dataRaw/` contains three sub-directories: `methylation/`, which contains files labeled `{sample_id}.bed`; `variant-calls/`, which contains VCF files similarly named `{sample_id}.pass.vcf`; and `phenotypes/` which contains a CSV with sample covariates such as age, sex, BMI, and disease status. `sample_id`s are three-digit numbers.

`dataDerived/` contains the processed (encoded) data. Files are named as `{sample_id}.{chromosome}.{encoding_scheme}.pt`. Briefly, these files contain tensors 

`scripts/` contains R and Python scripts to munge data. It would be best practice to make this all python, but I'm sticking with what I know for now.