# compute_cpg_summaty_stats.py
# For each chromosome, compute the summary statistics
# These stats will only be computed on the files found in the train/valid/test dirs
# In other words, IQR represents IQR for the 185 train samples (NOT all 381 samples)

import os
import argparse
import pickle
import pandas as pd

def load_as_df(ff):
    '''Helper function to load responses from pickle'''
    file = open(ff, 'rb')
    # dump information to that file
    data = pickle.load(file)
    
    # close the file
    file.close()
    tmp = pd.DataFrame.from_dict(
        data, orient = "index", 
        columns = ['methylated', 'coverage']
    )
    tmp['percentage'] = tmp['methylated'] / tmp['coverage']
    return(tmp[['percentage']])



def main(args):
    print(args)
    idir = args.idir 
    odir = args.odir

    all_substrings = [x.split('.')[1] for x in os.listdir(idir)]
    all_chroms = set([x for x in all_substrings if 'chr' in x])
    
    if not os.path.exists(odir):
        os.makedirs(odir)

    for chrom in all_chroms:
        print(f"Working on {chrom}")
        ofile = os.path.join(args.odir, f"{chrom}.stats.csv")
        all_dfs = []
    
        for ix, ff in enumerate(os.listdir(idir)):
            if ff.endswith(f".{chrom}.response.pkl"):
                df = load_as_df(os.path.join(idir, ff))
                df['subject'] = len(all_dfs) # arbitrary ID
                all_dfs.append(df)
                
                if len(all_dfs) % 25 == 0:
                    print(f"Loaded {len(all_dfs)} datasets so far")
    
        full_df = pd.concat(all_dfs)
    
        print("Pivoting to wide, KEEPING missing")
        wide = full_df.pivot(columns = 'subject', values = 'percentage')
    
        # Compute all the summary stats in vectorized manner
        out_df = pd.DataFrame(index = wide.index)
        
        print("Computing summary stats")
        print("   ...mean, min, max")
        out_df['mean'] = wide.mean(axis=1, skipna=True)
        out_df['min'] = wide.min(axis=1, skipna=True)
        out_df['max'] = wide.max(axis=1, skipna=True)
        out_df['std'] = wide.std(axis=1, skipna=True)
    
        print("   ...quantiles")
        out_df['q25'] = wide.quantile(q=0.25, axis=1)
        out_df['q50'] = wide.quantile(q=0.5, axis=1)
        out_df['q75'] = wide.quantile(q=0.75, axis=1)
        out_df['n_nonmissing'] = wide.count(axis=1)
        out_df['n_total'] = len(wide.columns)
    
        print("Applying final touches")
        out_df.reset_index(inplace=True)
        out_df = out_df.rename(columns = {'index':'position'})
        out_df['chrom'] = chrom
    
        out_df.to_csv(ofile, index = False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--idir', help = "Where the input data is stored")
    parser.add_argument('--odir', default = "dataDerived/statsForCpGs/")
    args = parser.parse_args()

    main(args)