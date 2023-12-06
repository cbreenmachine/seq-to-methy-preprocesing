import argparse
import pandas as pd
import os
import random

# Training and testing scheme
train_test_dict = {
    "group": ["train", "test"]*10 + ["valid"]*2,
    "chrom": ["chr" + str(x) for x in range(1, 23)]
}

train_test_df = pd.DataFrame(train_test_dict)

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

    # Create sample, group (100, 'train')
    sample_df = pd.DataFrame({'sample': samples[0:N], 'group': group})
    
    out_df = sample_df.merge(train_test_df, how= "outer", on="group")
    out_df.to_csv(args.ofile, index = False)
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--idir', default = "dataRaw/methylation/")
    parser.add_argument('--ofile', default = "dataDerived/randomization.csv")
    
    parser.add_argument('--n_train', help = "number of train samples", 
                        type=int, default = 10)
    parser.add_argument('--n_valid', help = "number of validation samples", 
                        type=int, default = 10)
    parser.add_argument('--n_test', help = "number of test samples", 
                        type=int, default = 10)

    args = parser.parse_args()
    main(args)
    