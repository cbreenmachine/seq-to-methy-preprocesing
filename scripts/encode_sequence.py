import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse
import torch

from encoder import OneHotEncoder, VariantEncoder, BivariantEncoder, AdditiveEncoder
from encoder import infer_chrom, infer_sample, load_vcf, load_bed, load_fa

# Changes
# 1. Read assignment from randomization.csv and write to train/, valid/, or test/


def main(args):

    # Define chrom/sample name
    chrom = infer_chrom(args.ofile_prefix)
    sample = infer_sample(args.ofile_prefix)
        
    # Read in data
    ref_seq = load_fa(f"{args.ref_dir}/{chrom}.fa", chrom)
    snp_df = load_vcf(f"{args.snp_dir}/{sample}.pass.vcf", chrom)
    methy_df = load_bed(f"{args.methy_dir}/{sample}.bed", chrom)

    # One-hot encoding first
    print("Encoding one-hot reference...")
    scheme1 = OneHotEncoder(ref_seq, methy_df, snp_df)
    scheme1.encode_data()
    out1 = torch.tensor(scheme1.encoded_tensor, dtype = torch.uint8)
    torch.save(out1, f"{args.ofile_prefix}.reference.pt")

    # Variant encoding
    print("Encoding warm variants...")
    scheme2 = VariantEncoder(ref_seq, methy_df, snp_df)
    scheme2.encode_data()
    out2 = torch.tensor(scheme2.encoded_tensor, dtype = torch.float)
    torch.save(out2, f"{args.ofile_prefix}.variant.pt")

    # Bi-variant (extended state space)
    print("Encoding one-hot bivariants...")
    scheme3 = BivariantEncoder(ref_seq, methy_df, snp_df)
    scheme3.encode_data()
    out3 = torch.tensor(scheme3.encoded_tensor, dtype = torch.uint8)
    torch.save(out3, f"{args.ofile_prefix}.bivariant.pt")

    # SNP mask (statistical genetics)
    print("Encoding number of alt alleles...")
    scheme4 = AdditiveEncoder(ref_seq, methy_df, snp_df)
    scheme4.encode_data()
    out4 = torch.tensor(scheme4.encoded_tensor, dtype = torch.uint8)
    torch.save(out4, f"{args.ofile_prefix}.additive.pt")


    # Keep track of positions, methylated counts, total counts as well
    # This file will be relatively small
    methy_df.to_pickle(f"{args.ofile_prefix}.responses.pkl")
    


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_dir', default = "dataRaw/reference/")
    parser.add_argument('--methy_dir', default = "dataRaw/methylation/")
    parser.add_argument('--snp_dir', default = "dataRaw/variant-calls/")
    parser.add_argument('--randomization', default = "dataRaw/randomization.csv")
    
    parser.add_argument('--sample', default = "101")
    args = parser.parse_args()
    
    main(args)