import argparse
import pandas as pd

def main(args):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--idir', default = "dataDerived/statsForCpGs/", help = "directory with files like chr17.stat.csv")
    parser.add_argument('--percent_nonmissing', default = 0.5)
    parser.add_argument('--min_iqr', default = 0.05)

    args = parser.parse_args()

    main(args)