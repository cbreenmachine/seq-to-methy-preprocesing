import pandas as pd
import argparse



def process_manifest(ifile):
     # Read CSV file into a DataFrame
    manifest = pd.read_csv(ifile, 
                           low_memory = False, skiprows=7, 
                           usecols=["IlmnID", "CHR", "MAPINFO"]).dropna()
    manifest['cpg_chromosome'] = 'chr' + manifest['CHR'].astype(str)
    manifest['cpg_position'] = manifest['MAPINFO']
    manifest = manifest.drop(columns = ["CHR", "MAPINFO"])

    return manifest


def process_meQTL(ifile):
    # meQTL data from 450k
    df = pd.read_csv(ifile, usecols=["cpg", "snp", "beta_a1", "allele1", "allele2", "freq_a1"])
    df['allele1'] = df['allele1'].str.upper()
    df['allele2'] = df['allele2'].str.upper()
    df['variant_call'] = "1/1"
    
    # Need coordinates split up a bit
    df[['snp_chromosome', 'snp_position', 'snp_type']] = df['snp'].str.split(':', n=2, expand = True)
    return df


def merge_and_filter(meQTL_df, manifest_df, max_dist = 1000):
    combined = meQTL_df.merge(manifest_df, left_on = 'cpg', right_on = 'IlmnID', how = 'inner').dropna()
    
    filtered = combined[combined['cpg_chromosome'] == combined['snp_chromosome']]
    filtered = filtered[filtered['snp_type'] == "SNP"]
    filtered['distance'] = filtered['cpg_position'].astype(int) - filtered['snp_position'].astype(int)
    
    final = filtered[filtered['distance'].abs() < max_dist]
    return final


def main(args):

    print("Loading manifest")
    manifest = process_manifest(args.manifest)

    print("Loading meQTLs")
    meqtl = process_meQTL(args.meqtl)

    print("Merging and filtering")
    final = merge_and_filter(meqtl, manifest)    
    
    # Write DataFrame to CSV
    final.to_csv(args.ofile, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--manifest', 
                        default = "data/meQTL/humanmethylation450_15017482_v1-2.csv",
                        help = "Manifest file provided by Illumina")
    parser.add_argument('--meqtl', default = "data/meQTL/assoc_meta_all.csv")
    parser.add_argument('--ofile', default = "data/meQTL/meQTL-variants.csv")

    args = parser.parse_args()

    main(args)


