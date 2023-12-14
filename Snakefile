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
# Reference genome URL
ref_genome_url = "ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

db_snp_url = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz"

FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

chrom_range = ["chr" + str(x) for x in range(1, 23)]
chrom_range_str = ",".join(chrom_range)

encode_samples = ['ENCFF534YXW', 'ENCFF933OKK', 'ENCFF573LIB', 'ENCFF844RZP', 'ENCFF068PZH']

###########################################################################
#################### Randomize to train and test ##########################
###########################################################################

df = pd.read_csv("dataDerived/randomization.csv")
samples_by_group = df.groupby('group')['sample'].apply(list).to_dict()

sample_range = df['sample'].tolist()


def get_ofiles_list_by_method(output_type, group, ext = "pt", odir="dataDerived/"):
    '''Handles the nested structure of the train/test samples'''
    # method: variant, reference
    samples = samples_by_group[group]
    chroms = chroms_by_group[group]

    files = [f"{s}.{c}.{output_type}.{ext}" for s,c in product(samples, chroms)]
    out = [os.path.join(odir, group, z) for z in files]
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
        get_ofiles_list_by_method('response', 'train', 'pkl'),
        get_ofiles_list_by_method('reference', 'train'),
        get_ofiles_list_by_method('variant', 'train'),
        get_ofiles_list_by_method('mask', 'train'),

        get_ofiles_list_by_method('response', 'valid', 'pkl'),
        get_ofiles_list_by_method('reference', 'valid'),
        get_ofiles_list_by_method('variant', 'valid'),
        "dataDerived/vcf.stats.csv",
        "dataDerived/vcf.stats.encode.csv",
        "dataDerived/covariates.csv"
       


###########################################################################
########################## VCF Summaries ##################################
###########################################################################


rule download_common_snps:
    input: HTTP.remote(db_snp_url, keep_local = False)
    output: "dataRaw/dbSNP/common_all_20180418.vcf.gz"
    conda: "envs/bcftools.yaml"
    shell: 
        """
        mv {input} {output}
        bcftools index {output}
        """
        

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
        chrom_range_str
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
    params:
        chrom_range_str
    conda: "envs/bcftools.yaml"
    shell: 
        """
        bcftools stats --regions {params[0]} {input[0]} {input[2]} > {output}
        """



rule combine_vcf_overlap_stats:
    input: 
        expand("dataDerived/variantComparisons/{sample}.out", sample = sample_range)
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

rule download_reference_genome:
    input:
        FTP.remote(ref_genome_url, keep_local=False)
    output:
        "dataRaw/reference/hg38.fna"
    shell:
        "gunzip -c {input} > {output}"

rule split_reference_by_chrom:
    input:
        "dataRaw/reference/hg38.fna"
    output:
        expand("dataRaw/reference/{chrom}.fa", chrom = chrom_range)
    shell:
        "python scripts/split_ref_seq_by_chrom.py --ifile {input}"

    

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

###########################################################################
############## Tell me when jobs are done or break ########################
###########################################################################

onsuccess: shell("mail -s 'DONE' {email_address} < {log}")
onerror: shell("mail -s 'Error for sequence processing pipeline' {email_address} < {log}")
