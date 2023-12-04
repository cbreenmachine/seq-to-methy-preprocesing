from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider



email_address = "cebreen@wisc.edu"

###########################################################################
############################## Constants ##################################
###########################################################################
# Reference genome URL
ref_genome_url = "ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

# 110 total
n_train = 55
n_valid = 5
n_test = 50

FTP = FTPRemoteProvider()

# Placeholder for now
chrom_range = ['chr' + str(x) for x in range(16, 23)]

# The three output types produced by encode_...py
suffixes = ['reference', 'variant', 'bivariant', 'additive']


###########################################################################
############################## Last file ##################################
###########################################################################
rule all:
    input: 
        "dataDerived/randomization.csv",
        expand("dataRaw/reference/{chrom}.fa", chrom = chrom_range)
    # input: 
    #     expand("data/cleaned/{sample}.{chrom}.{out}.pt", sample = sample_range, chrom = chrom_range, out = suffixes)

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


rule randomly_select_samples:
    params:
        n_train = n_train, 
        n_valid = n_valid, 
        n_test = n_test
    output: "dataDerived/randomization.csv"
    shell:
        """
        python scripts/randomize_samples.py \
        --n_train {params.n_train} \
        --n_valid {params.n_valid} \
        --n_test {params.n_test} \
        --ofile {output}
        """
    

# rule encode_samples:
#     params:
#         "data/cleaned/{sample}.{chrom}"
#     output:
#         "data/cleaned/{sample}.{chrom}.reference.pt",
#         "data/cleaned/{sample}.{chrom}.variant.pt",
#         "data/cleaned/{sample}.{chrom}.bivariant.pt",
#         "data/cleaned/{sample}.{chrom}.additive.pt"
#     resources:
#         mem_mb = 20000
#     shell:
#         "python scripts/encode_sequence.py --ofile_prefix {params}"

###########################################################################
############## Tell me when jobs are done or break ########################
###########################################################################
onsuccess:
   shell("mail -s 'DONE' {email_address} < {log}")

onerror:
   shell("mail -s 'Error for sequence processing pipeline' {email_address} < {log}")

    
        
        
        
        
