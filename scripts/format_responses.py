import pandas as pd
import argparse
import torch
import pickle
from utils.read_write import load_methylation_bed
from utils.constants import chroms_by_group

def make_input_name(idir, sample):
    return f"{idir}/{sample}.bed"

def make_output_name(odir, sample, chrom):
    return f"{odir}/{sample}.{chrom}.response.pkl"    


def format_to_dict(df):
    df = df.copy()
    df['response'] = df.loc[:, ['methylated', 'coverage']].apply(tuple, axis=1)
    return df['response'].to_dict()

    
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
        # Chromosome specific
        ofile_name = make_output_name(odir, sample, chrom)
        
        methy_chrom_df = methy_df.loc[methy_df['chrom'] == chrom, :]
        data = format_to_dict(methy_chrom_df)

        pickle_file = open(ofile_name, 'wb')
        pickle.dump(data, pickle_file)
        pickle_file.close()

        print(f"Wrote {ofile_name}")
    


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--idir', default = "dataRaw/methylation/")
    parser.add_argument('--odir', default = "dataDerived/")
    parser.add_argument('--group', default = "train", choices = ['train', 'valid', 'test'])
    parser.add_argument('--sample', default = "101")
    
    args = parser.parse_args()
    
    main(args)