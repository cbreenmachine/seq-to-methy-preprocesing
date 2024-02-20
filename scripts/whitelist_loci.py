import argparse
import pandas as pd
import os
import pickle

def main(args):

    output = {}

    for x in range(1, 23):
        ifile = os.path.join(args.idir, f"chr{x}.stats.csv")
        df = pd.read_csv(ifile)

        df['p_nonmissing'] = df['n_nonmissing'] / df['n_total']
        df['IQR'] = df['q75'] - df['q25']
        
        keepix = (df['p_nonmissing'] >= args.percent_nonmissing) & (df['IQR'] >= args.min_iqr)
        df_sub = df.loc[keepix]

        output[f"chr{x}"] = df_sub['position'].tolist()

    
    with open(args.ofile, 'wb') as handle:
        pickle.dump(output, handle, protocol=pickle.HIGHEST_PROTOCOL)    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--idir', default = "dataDerived/statsForCpGs/", help = "directory with files like chr17.stat.csv")
    parser.add_argument('--percent_nonmissing', default = 0.5)
    parser.add_argument('--min_iqr', default = 0.05)
    parser.add_argument('--ofile', default = "dataDerived/whitelists/IQR05.pkl")

    args = parser.parse_args()

    main(args)