import pandas as pd

ifile = "../data/phenotypes/2023-10-03-master-samplesheet.csv"
ofile = "../data/cleaned/covariates.csv"


def main(args):

    # Load data from CSV
    data = pd.read_csv(args.ifile)

    # Perform processing
    data = data.assign(
        group=lambda x: pd.np.where(x['diagnostic_group'] == "LOAD", 1,
                                    pd.np.where(x['diagnostic_group'] == "MCI", 0.5, 0)),
        bmi=lambda x: x['bmi'] / x['bmi'].max(skipna=True),
        age=lambda x: x['age_at_visit'] / x['age_at_visit'].max(skipna=True)
    )

    # Replacing one missing BMI value
    data['bmi'] = data['bmi'].fillna(0.5)

    # Write out the result
    data[['sample', 'group', 'bmi', 'age']].to_csv(args.ofile, index=False)


if __init__ == "__main__":
    