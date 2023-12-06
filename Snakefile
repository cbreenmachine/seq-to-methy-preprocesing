from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
import pandas as pd
from scripts.utils.read_write import get_encoding_file_names, chroms_by_method

email_address = "cebreen@wisc.edu"


###########################################################################
############################## Constants ##################################
###########################################################################
# Reference genome URL
ref_genome_url = "ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
FTP = FTPRemoteProvider()
chrom_range = ["chr" + str(x) for x in range(1, 23)]

###########################################################################
#################### Randomize to train and test ##########################
###########################################################################

df = pd.read_csv("dataDerived/randomization.csv")
samples_by_method = df.groupby('group')['sample'].apply(list).to_dict()


def generate_sparse_output_files(wildcards, samples_by_method, chroms_by_method):
    out = []

    for g in ['train', 'valid', 'test']:
        out.append(
            expand(
                "dataDerived/{group}/{sample}.{chrom}.variant.pt", 
                group = g,
                sample = samples_by_method[g],
                chrom = chroms_by_method[g]
            )
        )
    return(out)


def get_sample(wildcards, output):
    bn = os.path.basename(output[0])
    return bn.split(".")[0]

def get_chrom(wildcards, output):
    bn = os.path.basename(output[0])
    return bn.split(".")[1]

def get_group(wildcards, output):
    no_file = os.path.dirname(output[0])
    return no_file.split("/")[-1]
    
    


###########################################################################
############################## Last file ##################################
###########################################################################
rule all:
    input: 
        expand("dataDerived/train/{sample}.{chrom}.variant.pt", 
                sample = samples_by_method['train'], 
                chrom = chroms_by_method['train']),

        expand("dataDerived/valid/{sample}.{chrom}.variant.pt", 
                sample = samples_by_method['valid'], 
                chrom = chroms_by_method['valid'])
        
        
        
   
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


rule encode_variant_samples:
    params:
        sample = lambda wildcards, output: get_sample(wildcards, output),
        group = lambda wildcards, output: get_group(wildcards, output)
    output: 
        expand("dataDerived/{{group}}/{{sample}}.{chrom}.variant.pt", chrom = chrom_range)
    resources:
        mem_mb = 2000
    log: "logs/encode_variant_samples.log"
    benchmark: "benchmarks/encode_variant_samples.benchmark"
    shell:
        """
        python scripts/encode_sequence.py \
            --sample {params.sample} \
            --group {params.group} \
            --method variant
        """

###########################################################################
############## Tell me when jobs are done or break ########################
###########################################################################

onsuccess: shell("mail -s 'DONE' {email_address} < {log}")
onerror: shell("mail -s 'Error for sequence processing pipeline' {email_address} < {log}")

    
        
        
        
        
