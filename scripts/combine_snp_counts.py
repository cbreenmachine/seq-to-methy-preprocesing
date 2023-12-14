from os.path import basename
from os import listdir
import argparse
import pandas as pd


def extract_counts(input_file):
    '''Parse bcftools stats output'''
    out = []
    out.append(basename(input_file).split(".")[0])
    
    with open(input_file, 'r') as file:
        for line_number, line in enumerate(file, 1):
            if "number of SNPs:" in line:
                words = line.split()
                count = words[5]
                out.append(count)
    return out


def main(args):

    all_files = listdir(args.idir)
    data = []

    for input_file in all_files:
        if ".out" in input_file:
            data.append(extract_counts(f"{args.idir}/{input_file}"))

    df = pd.DataFrame(data, columns = ['sample', 'sample_only', 'dbSNP_only', 'intersection'])
    df.to_csv(args.ofile, index = False)




if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--idir', default = "dataDerived/variantComparisons/vcf.stats.csv")
    parser.add_argument('--ofile', default = "dataDerived/vcf.stat.csv")
    args = parser.parse_args()
    
    main(args)