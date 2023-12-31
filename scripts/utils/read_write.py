import numpy as np
import pandas as pd
from Bio import SeqIO



def get_encoding_file_names(rando_file_path, extension, group, odir = "dataDerived"):
    df = pd.read_csv(rando_file_path)

    out = []
    for _, row in df.iterrows():
        g, s = row['group'], row['sample']

        if g == group:
            for c in chroms_by_method[g]:
                out.append(f"{odir}/{g}/{s}.{c}.{extension}.pt")
        
    return(out)


def get_encoding_file_prefix(rando_file_path, extension, group, odir = "dataDerived"):
    df = pd.read_csv(rando_file_path)

    out = []
    for _, row in df.iterrows():
        g, s = row['group'], row['sample']

        if g == group:
            for c in chroms_by_method[g]:
                out.append(f"{odir}/{g}/{s}.{c}.{extension}.pt")
        
    return(out)


def infer_sample(x):
    '''Return the sample we're operating on'''
    return(x.split("/")[-1].split(".")[0])


def infer_chrom(x):
    '''Return the chromosome we're operating on'''
    return(x.split("/")[-1].split(".")[1])


def load_fa(ifile, chrom):
    '''wrapper to read from fasta'''
    return(SeqIO.to_dict(SeqIO.parse(open(ifile), 'fasta'))[chrom])



def load_vcf(ifile):
    '''
    Reads and extracts relevant information from VCF file

    Returns
    pd.DataFrame
    chrom pos reference alternate variant_call
    chr1  23  A         G         0/0
    '''
    df = pd.read_csv(ifile, 
                     sep='\t', 
                     comment='#', 
                     header=None, 
                     usecols=[0, 1, 3, 4, 9],
                     names=['chrom', 'pos', 'reference', 'alternate', 'extra_info']
                     )

    return(df)


def subset_vcf(df, chrom):
    
    df = df.copy()
    
    # Filter to shrink the number of comparisons that need to be made
    df = df.loc[df['chrom'] == chrom]
    
    # Pull out the SNP call
    df['variant_call'] = df['extra_info'].str[:3]

    # Get rid of Ns, indicate that ref homozygous
    df = df[df['reference'].isin(['A', 'C', 'G', 'T'])]
    
    df['data'] = df[['reference', 'alternate', 'variant_call']].apply(tuple, axis=1)
    df = df.set_index("pos")
    
    return df['data'].to_dict()

def load_methylation_bed(ifile, min_coverage = 10):
    '''
    
    '''
    full_df = pd.read_table(ifile)

    keepix = ((full_df['coverage'] >= min_coverage) &
              (full_df['methylated'] <= full_df['coverage']))
    
    sub_df = full_df.loc[keepix].copy(deep=True)

    # Convert 0-based BED coordinate to 1-based VCF type
    sub_df['pos'] = sub_df.chromStart + 1
    sub_df.set_index('pos', inplace = True)
    
    return(sub_df[['chrom', 'methylated', 'coverage', 'strand']])

