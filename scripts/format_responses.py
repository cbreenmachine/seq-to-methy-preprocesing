import pandas as pd
import argparse
import torch
from utils.read_write import load_methylation_bed, chroms_by_group

def make_input_name(idir, sample):
    return f"{idir}/{sample}.bed"

def make_output_name(odir, sample, chrom):
    return f"{odir}/{sample}.{chrom}.responses.pkl"    

    
def main(args):

    # Define chrom/sample name
    sample = args.sample
    group = args.group
    idir = args.idir

    # Output directory incorporates train/test
    odir = f"{args.odir}/{args.group}/"

    # Load load load
    methy_input_file = make_input_name(idir, sample)    
    methy_df = load_methylation_bed(methy_input_file)
    
    for chrom in chroms_by_group[group]:
        methy_chrom_df = methy_df[methy_df['chrom'] == chrom]

        # Chromosome specific
        output_file_name = make_output_name(odir, sample, chrom)
        methy_chrom_df.to_pickle(output_file_name)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--idir', default = "../dataRaw/methylation/")
    parser.add_argument('--odir', default = "../dataDerived/")
    parser.add_argument('--group', choices = ['train', 'valid', 'test'])
    parser.add_argument('--sample', default = "101")
    
    args = parser.parse_args()
    
    main(args)