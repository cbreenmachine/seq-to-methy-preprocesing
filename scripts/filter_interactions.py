#!/usr/bin/env python

import argparse
from liftover import get_lifter
import pandas as pd
import os

# CONSTANTS
type_dict = {'baitChr': 'str', 'oeChr': 'str'}
chrom_range = [str(x) for x in range(1, 23)]
col_names = ['chrom', 'start', 'end', 'interaction_id']


converter = get_lifter('hg19', 'hg38')

def main(args):

    promoter_intervals = []
    enhancer_intervals = []

    inter_df = pd.read_table(args.ifile, dtype = type_dict)
    
    for chrom in chrom_range:
        
        sub_df = inter_df.loc[inter_df.baitChr == chrom]
        
        for ix, row in sub_df.iterrows():
    
            # Pluck out baits (proxy for promoters)
            bait_chrom = row['baitChr']
            oe_chrom = row['oeChr']
    
            if bait_chrom == oe_chrom:
    
                # Pull out the start and end for baits and enhancers repsectively
                bait_start = row['baitStart']
                bait_end = row['baitEnd']
    
                oe_start = row['oeStart']
                oe_end = row['oeEnd']
            
                # Convert start and end coordinate separately
                bait_left = converter[bait_chrom][bait_start]
                bait_right = converter[bait_chrom][bait_end]
            
                # Convert start and end coordinate separately
                oe_left = converter[oe_chrom][oe_start]
                oe_right = converter[oe_chrom][oe_end]
    
                # Verify that all conversions worked
                if (
                    len(bait_left) > 0 and len(bait_right) > 0 and 
                    len(oe_left) > 0 and len(oe_right) > 0
                ):
                    # When we want to match these enhancers/promoters up later
                    interaction_id = chrom + "-" + str(ix)
                    promoter_intervals.append(
                        (bait_chrom, bait_left[0][1], bait_right[0][1], interaction_id)
                    )
                    enhancer_intervals.append(
                        (oe_chrom, oe_left[0][1], oe_right[0][1], interaction_id)
                    )
    
        # Convert to data frame for writing csv
        promoters = pd.DataFrame(promoter_intervals, columns = col_names)
        enhancers = pd.DataFrame(enhancer_intervals, columns = col_names)
        
        # Create names
        ofile_p = os.path.join(args.odir, "promoters.chr{0}.csv".format(chrom))
        ofile_e = os.path.join(args.odir, "enhancers.chr{0}.csv".format(chrom))
        
        promoters.to_csv(ofile_p, index=False)
        enhancers.to_csv(ofile_e, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='filter_interactions',
        description='filters and lifts interactions'
    )
    
    parser.add_argument('--ifile',  default = "../data/interactions.txt")
    parser.add_argument('--odir',  default = "../data/")
    
    args = parser.parse_args()
    main(args)