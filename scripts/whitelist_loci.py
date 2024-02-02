import argparse
import pandas as pd

def main(args):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--', 
                        default = "data/meQTL/humanmethylation450_15017482_v1-2.csv",
                        help = "Manifest file provided by Illumina")
    parser.add_argument('--meqtl', default = "data/meQTL/assoc_meta_all.csv")
    parser.add_argument('--ofile', default = "data/meQTL/meQTL-variants.csv")

    args = parser.parse_args()

    main(args)