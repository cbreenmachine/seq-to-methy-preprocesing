import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse
import torch

from vencoder.onehot import OneHotEncoder
from vencoder.variant import VariantEncoder

from utils.read_write import load_vcf, load_bed, load_fa, get_encoding_file_names, chroms_by_method

def encoder_factory(method):
    if method == "variant":
        return VariantEncoder
    elif method == "reference":
        return OneHotEncoder


def main(args):

    # Define chrom/sample name
    sample = args.sample
    method = args.method
    group = args.group
    
    ref_seq_dir = args.ref_seq_dir
    var_calls_dir = args.var_calls_dir
    
    odir = f"{args.odir}/{args.group}/"

    for chrom in chroms_by_method[group]:
        # Derived some values
        ofile = f"{odir}/{sample}.{chrom}.{method}.pt"
        print(f"Will output {ofile}")
        
        # Read in sequence data and convert to one long string
        ref_sequence = load_fa(f"{ref_seq_dir}/{chrom}.fa", chrom)
        var_calls = load_vcf(f"{var_calls_dir}/{sample}.pass.vcf", chrom)

        # Use the factory method to refactor
        encoder_instance = encoder_factory(method)
        my_encoder = encoder_instance(ref_sequence, var_calls)
        my_encoder.encode_data()

        output_tensor = torch.tensor(my_encoder.encoding, dtype = torch.uint8)    
        torch.save(output_tensor, ofile)

    # Keep track of positions, methylated counts, total counts as well
    # This file will be relatively small

    # response_output_file = f"{args.ofile_prefix}.responses.pkl"
    # if not os.file.exists(response_output_file):
    #     methy_df = load_bed(f"{args.methy_dir}/{sample}.bed", chrom)
    #     methy_df.to_pickle(response_output_file)
    


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_seq_dir', default = "dataRaw/reference/")
    parser.add_argument('--var_calls_dir', default = "dataRaw/variant-calls/")
    parser.add_argument('--odir', default = "dataDerived/")

    parser.add_argument('--group', choices = ['train', 'valid', 'test'])
    parser.add_argument('--sample', default = "101")
    parser.add_argument('--method', choices = ["reference", "variant"])
    
    args = parser.parse_args()
    
    main(args)