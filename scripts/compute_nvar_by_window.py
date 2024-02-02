import sys
sys.path.append("..")
import pandas as pd
from itertools import product
import argparse

from snvme.datasets.munger import DataFileNames, MethyFileOpener
from snvme.datasets.generators import expressive_generator

# Randomly selected 10 samples
samples = ["319", "297", "389", "153", "123", "140", "457", "415", "493", "388"]
chroms = ["chr1"]

data = []

def main(args):

    for s, c in product(samples, chroms):
        for window in [500, 1000, 2000, 4000, 8000]:
            print(f"Window size: {window}")
            
            common = DataFileNames(
                root_dir = args.idir, group = "train",
                type = "reference", window = window
            )
        
            methy_pipe = MethyFileOpener(
                common, s, c, False, 
                sample_prob_no_variants = 1,
                sample_prob_with_variants = 1
            )
            
            methy_pipe.load_data()
        
            for  _, extra in expressive_generator(methy_pipe):
                data.append([window, extra['n_var'], s, c])

        df = pd.DataFrame(data, columns=['window_size', 'n_variants', 'sample', 'chrom'])
        
        try:
            df.to_csv(args.ofile, index=False)
        except OSError:
            df.to_csv("nvar_by_window.csv", index = False)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--idir', default = "../../seq-to-methy-preprocesing/dataDerived/")
    parser.add_argument('--ofile', default = "../data/variant_summaries/nvar_by_window.csv")
    args = parser.parse_args()

    main(args)


