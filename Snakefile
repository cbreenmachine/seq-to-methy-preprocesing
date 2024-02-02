from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import pandas as pd
from scripts.utils.constants import chroms_by_group

from itertools import product
import glob

email_address = "cebreen@wisc.edu"

###########################################################################
############################## Constants ##################################
###########################################################################

all_chroms = ["chr" + str(x) for x in range(1, 23)]
all_chroms_str = ",".join(all_chroms)

encode_samples = ['ENCFF534YXW', 'ENCFF933OKK', 'ENCFF573LIB', 'ENCFF844RZP', 'ENCFF068PZH']




###########################################################################
#################### Randomize to train and test ##########################
###########################################################################

df = pd.read_csv("dataDerived/randomization.csv")
all_samples = df['sample'].tolist()

# Samples by group
samples_by_group = df.groupby('group')['sample'].apply(list).to_dict()

train_samples = samples_by_group['train']
valid_samples = samples_by_group['valid']
test_samples = samples_by_group['test']

# Chroms by group
train_chroms = ["chr" + str(x) for x in range(1, 20, 2)]
valid_chroms = ["chr21", "chr22"]
test_chroms = ["chr" + str(x) for x in range(2, 21, 2)]

def get_ofiles_list_by_method(samples, chroms, odir):
    '''Handles the nested structure of the train/test samples

    Outputs like sample.chr.variant.pt
    '''

    out = []

    for ext in ['variant.pt', 'reference.pt', 'mask.pt', 'response.pkl']:
        for s, c in product(samples, chroms):
            out.append(os.path.join(odir, f"{s}.{c}.{ext}"))

    return(out)



###########################################################################
############################## Last file ##################################
###########################################################################


def make_encoding_input_names(wildcards):
    out = []

    for m in ['reference', 'variant', 'mask']:
        for g in ['train', 'valid']:
            out += get_ofiles_list_by_method(m, g)
    return(out)
    


rule all:
    input: 
        get_ofiles_list_by_method(train_samples, train_chroms, "dataDerived/train/"),
        get_ofiles_list_by_method(valid_samples, valid_chroms, "dataDerived/valid/"),
        get_ofiles_list_by_method(test_samples, test_chroms, "dataDerived/test/"),

        "dataDerived/vcf.stats.csv",
        "dataDerived/vcf.stats.encode.csv",
        "dataDerived/covariates.csv",
        expand("dataDerived/statsForCpGs/{chr}.stats.csv", chr = train_chroms)
        



#######################################################################
###################### Public data downloads ##########################
#######################################################################

rule hg19:
    output: "dataPublic/hg19/hg19.fna"
    shell: "bash scripts/process)hg19.sh"

rule hg38:
    output: "dataPublic/hg38/split/chrY.fa"
    shell: "bash scripts/process_hg38.sh"

rule dbSNP:
    output: 
        "dataPublic/dbSNP/common_all_20180418.renamed.vcf.gz.csi"
    shell: "bash scripts/process_dbSNP.sh"

rule cpg_stats_train:
    output: 
        expand("dataDerived/statsForCpGs/{chrom}.stats.csv", chrom = train_chroms)
    shell:
        """
        python scripts/compute_cpg_summary_stats.py --idir dataDerived/train/
        """

rule cpg_stats_valid:
    output: 
        expand("dataDerived/statsForCpGs/{chrom}.stats.csv", chrom = valid_chroms)
    shell:
        """
        python scripts/compute_cpg_summary_stats.py --idir dataDerived/valid/
        """

rule cpg_stats_test:
    output: 
        expand("dataDerived/statsForCpGs/{chrom}.stats.csv", chrom = test_chroms)
    shell:
        """
        python scripts/compute_cpg_summary_stats.py --idir dataDerived/test/
        """

###########################################################################
########################## VCF Summaries ##################################
###########################################################################


        
rule process_common_snps:
    input:
        dbsnp_file = "dataRaw/dbSNP/common_all_20180418.vcf.gz",
        chr_map = "dataRaw/dbSNP/N_to_chrN.txt"
    output:
        "dataRaw/dbSNP/common_all_20180418.renamed.vcf.gz",
        "dataRaw/dbSNP/common_all_20180418.renamed.vcf.gz.csi"
    conda: "envs/bcftools.yaml"
    shell:
        """
        bcftools annotate --rename-chrs {input.chr_map} {input.dbsnp_file} | bgzip > {output[0]}
        bcftools index {output[0]}
        """


rule index_vcfs:
    input: "dataRaw/variants/{sample}.pass.vcf.gz"
    output: "dataRaw/variants/{sample}.pass.vcf.gz.csi"
    conda: "envs/bcftools.yaml"
    shell: "bcftools index -f {input} > {output}"



rule compress_and_index_encode_samples:
    conda: "envs/bcftools.yaml"
    input:
        "dataRaw/ENCODE/{sample}.vcf.gz"
    output:
        "dataRaw/ENCODE/{sample}.bcf",
        "dataRaw/ENCODE/{sample}.bcf.csi"
    shell:
        """
        bcftools view -Ob -o {output[0]} {input[0]}
        bcftools index {output[0]} 
        """


rule compute_vcf_overlaps:
    input: 
        "dataRaw/variants/{sample}.pass.vcf.gz",
        "dataRaw/variants/{sample}.pass.vcf.gz.csi",
        "dataRaw/dbSNP/common_all_20180418.renamed.vcf.gz",
        "dataRaw/dbSNP/common_all_20180418.renamed.vcf.gz.csi"
    output:
        "dataDerived/variantComparisons/{sample}.out"
    params:
        all_chroms_str
    conda: "envs/bcftools.yaml"
    shell: 
        """
        bcftools stats --regions {params[0]} {input[0]} {input[2]} > {output}
        """

rule compute_vcf_overlaps_encode:
    input: 
        "dataRaw/ENCODE/{sample}.bcf",
        "dataRaw/ENCODE/{sample}.bcf.csi",
        "dataRaw/dbSNP/common_all_20180418.renamed.vcf.gz",
        "dataRaw/dbSNP/common_all_20180418.renamed.vcf.gz.csi"
    output:
        "dataDerived/variantComparisonsENCODE/{sample}.out"
    params: all_chroms_str
    conda: "envs/bcftools.yaml"
    shell: 
        """
        bcftools stats --regions {params[0]} {input[0]} {input[2]} > {output}
        """


rule combine_vcf_overlap_stats:
    input: 
        expand("dataDerived/variantComparisons/{sample}.out", sample = all_samples)
    params: 
        "dataDerived/variantComparisons/"
    output:
        "dataDerived/vcf.stats.csv"
    shell:
        """
        python scripts/combine_snp_counts.py \
            --idir {params[0]} \
            --ofile {output}
        """



rule combine_vcf_overlap_stats_encode:
    input: 
        expand("dataDerived/variantComparisonsENCODE/{sample}.out", sample = encode_samples)
    params: 
        "dataDerived/variantComparisonsENCODE/"
    output:
        "dataDerived/vcf.stats.encode.csv"
    shell:
        """
        python scripts/combine_snp_counts.py \
            --idir {params[0]} \
            --ofile {output}
        """
    
    
   
###########################################################################
############################## Workflow ###################################
###########################################################################



rule clean_covariates:
    input: "dataRaw/phenotypes/2023-10-03-master-samplesheet.csv"
    output: "dataDerived/covariates.csv"
    shell: 
        """
        python scripts/standardize_covariates.py \
            --ifile {input} \
            --ofile {output}
        """


###########################################################################
################ Encode train, validation and test ########################
###########################################################################


rule encode_train_samples:
    priority: 10
    output: 
        expand("dataDerived/train/{{sample}}.{chr}.{{method}}.pt", 
                chr = chroms_by_group['train'])
    shell:
        """
        python scripts/encode_sequence.py \
            --sample {wildcards.sample} \
            --group train \
            --method {wildcards.method}
        """

rule encode_valid_samples:
    priority: 10
    output: 
        expand("dataDerived/valid/{{sample}}.{chr}.{{method}}.pt", 
                chr = chroms_by_group['valid'])
    shell:
        """
        python scripts/encode_sequence.py \
            --sample {wildcards.sample} \
            --group valid \
            --method {wildcards.method}
        """

rule encode_test_samples:
    priority: 10
    output: 
        expand("dataDerived/test/{{sample}}.{chr}.{{method}}.pt", 
                chr = chroms_by_group['test'])
    shell:
        """
        python scripts/encode_sequence.py \
            --sample {wildcards.sample} \
            --group test \
            --method {wildcards.method}
        """


rule format_train_responses:
    output: 
        expand("dataDerived/train/{{sample}}.{chr}.response.pkl",
                chr = chroms_by_group['train'])
    shell:
        """
        python scripts/format_responses.py \
            --sample {wildcards.sample} \
            --group train
        """

rule format_valid_responses:
    output: 
        expand("dataDerived/valid/{{sample}}.{chr}.response.pkl",
                chr = chroms_by_group['valid'])
    shell:
        """
        python scripts/format_responses.py \
            --sample {wildcards.sample} \
            --group valid
        """

rule format_test_responses:
    output: 
        expand("dataDerived/test/{{sample}}.{chr}.response.pkl",
                chr = chroms_by_group['test'])
    shell:
        """
        python scripts/format_responses.py \
            --sample {wildcards.sample} \
            --group test
        """

###########################################################################
############## Tell me when jobs are done or break ########################
###########################################################################

onsuccess: shell("mail -s 'DONE' {email_address} < {log}")
onerror: shell("mail -s 'Error for sequence processing pipeline' {email_address} < {log}")