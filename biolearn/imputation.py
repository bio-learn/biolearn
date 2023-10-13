import pandas as pd
from biolearn.util import get_data_file


def impute_from_standard(dnam, cpg_averages):
    """
    Impute all missing values using the averages from an external dataset

    :param dnam: DataFrame with samples as columns and cpg sites as rows.
    :param cpg_averages: Series containing reference averages for cpg sites.
    :return: DataFrame with missing values filled.
    """

    df_filled = dnam.apply(lambda col: col.fillna(cpg_averages))
    return df_filled


def impute_from_average(dnam):
    """
    Impute all missing values using the average from the current dataset

    :param dnam: DataFrame with samples as columns and cpg sites as rows.
    :return: DataFrame with missing values filled.
    """

    row_averages = dnam.mean(axis=1)
    df_filled = dnam.T.fillna(row_averages).T
    return df_filled


def hybrid_impute(dnam, cpg_source, required_cpgs, threshold=0.8):
    """
    Imputes missing values based on a threshold. For sites that have less data than
    the threshold requires all values for that site are replaced from source.
    Otherwise the average of existing values in the input data are used.

    :param dnam: DataFrame with samples as columns and cpg sites as rows.
    :param cpg_source: Series containing reference averages for cpg sites.
    :param required_cpgs: List of cpgs that need to be in the final dataset.
    :param threshold: Threshold for determining imputation strategy. Default is 0.8.
    :return: DataFrame with missing values filled.
    """
    # Drop rows below the threshold, these will be replaced entirely from cpg_source
    cpgs_below_threshold = dnam.notna().mean(axis=1) < threshold
    dnam = dnam.drop(dnam[cpgs_below_threshold].index)

    # Impute remaining rows using impute_from_average
    df_filled = impute_from_average(dnam)

    missing_cpgs_from_dataset = set(required_cpgs) - set(df_filled.index)
    missing_cpgs_from_source = [
        cpg for cpg in missing_cpgs_from_dataset if cpg not in cpg_source
    ]

    if missing_cpgs_from_source:
        raise ValueError(
            f"Tried to fill the following cpgs but they were missing from cpg_source: {missing_cpgs_from_source}"
        )

    for cpg in missing_cpgs_from_dataset:
        df_filled.loc[cpg] = cpg_source.loc[cpg]

    return df_filled.sort_index()


def biolearn_impute(dnam):
    biolearn_averages_file = get_data_file("biolearn_averages_450k.csv")
    df = pd.read_csv(biolearn_averages_file, index_col=0)
    biolearn_averages = df["average"]
    return impute_from_standard(dnam, biolearn_averages)
