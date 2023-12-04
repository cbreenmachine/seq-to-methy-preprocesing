from Bio import SeqIO
import os
import argparse

def main(input_file):
    ifile = args.ifile

    # Read in reference fasta and construct output directory name
    reference = SeqIO.parse(ifile, format = "fasta")
    output_dir = os.path.dirname(ifile)

    # Loop thru each contig (chromosomes + some chr1_random contigs)
    for contig in reference:
        output_file = os.path.join(output_dir, f"{contig.id}.fa")
        SeqIO.write(contig, output_file, format = "fasta")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--ifile', default = "../dataRaw/reference/hg38.fna")
    args = parser.parse_args()
    
    main(args)