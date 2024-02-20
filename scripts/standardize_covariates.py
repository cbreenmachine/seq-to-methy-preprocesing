# standardize_covariates.py
# Coleman Breen | cebreen@wisc.edu | Updated February 20, 2024
# Subsets phenotype data and standardizes covairates to be supported on 
# [0, 1] (in the case of continuous/numeric variables) and {-1, 0, 1}
# in the case of ordinal. Updated to remove duplicate rows.

import pandas as pd
import numpy as np
import argparse


def main(args):

    # Load data from CSV
    data = pd.read_csv(args.ifile)

    # Perform processing
    data = data.assign(
        group=lambda x: np.where(x['diagnostic_group'] == "LOAD", 1,
                                    np.where(x['diagnostic_group'] == "MCI", 0, -1)),
        is_male=lambda x: np.where(x['sex'] == "Male", 1, 0),
        bmi=lambda x: x['bmi'] / x['bmi'].max(skipna=True),
        age=lambda x: x['age_at_visit'] / x['age_at_visit'].max(skipna=True)
    )

    # Replacing one missing BMI value
    data['bmi'] = data['bmi'].fillna(0.5)
    data['sample'] = data['sample_id']

    data = data.groupby('sample').nth(0)

    # Write out the result
    data[['sample', 'group', 'is_male', 'bmi', 'age']].to_csv(args.ofile, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ifile', default = "dataRaw/phenotypes/2023-10-03-master-samplesheet.csv")
    parser.add_argument('--ofile', default = "dataDerived/covariates.csv")
    args = parser.parse_args()
    
    main(args)