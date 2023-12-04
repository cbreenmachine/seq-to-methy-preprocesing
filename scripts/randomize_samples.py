import argparse
import pandas as pd
import os
import random

def main(args):

    random.seed(919) # Raleigh
    
    # List all files (there will be ., ..)
    all_files = os.listdir(args.idir)

    # Lop off the extension and extra info
    all_prefixes = [x.split(".")[0] for x in all_files]

    # 1 gets rid of the empty string
    samples = sorted(set(all_prefixes))[1: ] 

    # Samples will now be shuffled
    random.shuffle(samples)

    # Now we construct a list of group assignments
    group = (["train"] * args.n_train) + \
    (["valid"] * args.n_valid) + \
    (["test"] * args.n_test)

    # Subset the samples and write out to a csv
    N = args.n_train + args.n_valid + args.n_test
    output = pd.DataFrame({'sample': samples[0:N], 'group': group})
    output.to_csv(args.ofile, index = False)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--idir', default = "dataRaw/methylation/")
    parser.add_argument('--ofile', default = "dataDerived/randomization.csv")
    
    parser.add_argument('--n_train', help = "number of train samples", type=int)
    parser.add_argument('--n_valid', help = "number of validation samples", type=int)
    parser.add_argument('--n_test', help = "number of test samples", type=int)

    args = parser.parse_args()
    main(args)
    