import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse
import torch

from vencoder.onehot import OneHotEncoder
from vencoder.variant import VariantEncoder
from vencoder.additive import AdditiveEncoder

from utils.read_write import load_vcf, subset_vcf, load_fa, get_encoding_file_names
from utils.constants import chroms_by_group

def encoder_factory(method):
    if method == "variant":
        return VariantEncoder
    elif method == "reference":
        return OneHotEncoder
    elif method == "mask":
        return AdditiveEncoder


def main(args):

    # Define chrom/sample name
    sample = args.sample
    method = args.method
    group = args.group
    
    ref_seq_dir = args.ref_seq_dir
    var_calls_dir = args.var_calls_dir
    
    odir = f"{args.odir}/{args.group}/"

    try:
        var_calls_df = load_vcf(f"{var_calls_dir}/{sample}.pass.vcf")
    except OSError:
        var_calls_df = load_vcf(f"{var_calls_dir}/{sample}.pass.vcf.gz")


    for chrom in chroms_by_group[group]:
        # Derived some values
        ofile = f"{odir}/{sample}.{chrom}.{method}.pt"
        print(f"Will output {ofile}")

        # Subset dataframe
        calls_by_chrom = subset_vcf(var_calls_df, chrom)
        
        # Read in sequence data and convert to one long string
        ref_sequence = load_fa(f"{ref_seq_dir}/{chrom}.fa", chrom)

        # Use the factory method to refactor
        encoder_instance = encoder_factory(method)
        my_encoder = encoder_instance(ref_sequence, calls_by_chrom)
        my_encoder.encode_data()

        # Storage saving
        if method == "variant":
            output_tensor = torch.tensor(my_encoder.encoding, dtype = torch.float)    
        else:
            output_tensor = torch.tensor(my_encoder.encoding, dtype = torch.uint8)
        
        torch.save(output_tensor, ofile)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_seq_dir', default = "dataRaw/reference/")
    parser.add_argument('--var_calls_dir', default = "dataRaw/variants/")
    parser.add_argument('--odir', default = "dataDerived/")

    parser.add_argument('--group', choices = ['train', 'valid', 'test'])
    parser.add_argument('--sample', default = "101")
    parser.add_argument('--method', choices = ["reference", "variant", "mask"])
    
    args = parser.parse_args()
    
    main(args)