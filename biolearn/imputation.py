import pandas as pd
from biolearn.util import get_data_file


def impute_from_standard(dnam, cpg_averages, cpgs_to_impute=None):
    """
    Impute all missing values using the averages from an external dataset

    :param dnam: DataFrame with samples as columns and cpg sites as rows.
    :param cpg_averages: Series containing reference averages for cpg sites.
    :param cpgs_to_impute: List of CpG sites to impute.
    :return: DataFrame with missing values filled.
    """
    if cpgs_to_impute:
        impute_rows = dnam.loc[cpgs_to_impute]
        impute_rows = impute_rows.apply(lambda col: col.fillna(cpg_averages))
        df_filled = dnam.combine_first(impute_rows)
    else:
        df_filled = dnam.apply(lambda col: col.fillna(cpg_averages))
    return df_filled


def impute_from_average(dnam, cpgs_to_impute=None):
    """
    Impute all missing values using the average from the current dataset

    :param dnam: DataFrame with samples as columns and cpg sites as rows.
    :param cpgs_to_impute: List of CpG sites to impute.
    :return: DataFrame with missing values filled.
    """

    # Protect input from mutation
    dnam_copy = dnam.copy()
    means = dnam_copy.mean(axis=1)

    if cpgs_to_impute:
        mask = dnam_copy.loc[cpgs_to_impute].isna()
        dnam_copy.loc[cpgs_to_impute] = dnam_copy.loc[cpgs_to_impute].where(
            ~mask, means[cpgs_to_impute], axis=0
        )
    else:
        dnam_copy = dnam_copy.where(dnam_copy.notna(), means, axis=0)

    return dnam_copy


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
