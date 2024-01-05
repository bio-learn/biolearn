import pandas as pd
from biolearn.util import get_data_file


def impute_from_standard(dnam, cpg_averages, cpgs_to_impute=None):
    """
    Impute all missing values in a DNA methylation dataset using the averages from an external dataset.

    Args:
        dnam (pd.DataFrame): DataFrame with samples as columns and CpG sites as rows.
        cpg_averages (pd.Series): Series containing reference averages for CpG sites.
        cpgs_to_impute (list of str, optional): List of CpG sites to impute. Missing cpgs will only be imputed if in this list.

    Returns:
        pd.DataFrame: DataFrame with missing values filled.
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
    Impute all missing values in a DNA methylation dataset using the average from the dataset itself.

    Args:
        dnam (pd.DataFrame): DataFrame with samples as columns and CpG sites as rows.
        cpgs_to_impute (list of str, optional): List of CpG sites to impute.

    Returns:
        pd.DataFrame: DataFrame with missing values filled.
    """
    # Protect input from mutation
    dnam_copy = dnam.copy()
    means = dnam_copy.mean(axis=1)

    if cpgs_to_impute:
        # Filter out non-existent CpG sites
        existing_cpgs = [
            cpg for cpg in cpgs_to_impute if cpg in dnam_copy.index
        ]

        # Apply imputation only to existing CpG sites
        mask = dnam_copy.loc[existing_cpgs].isna()
        dnam_copy.loc[existing_cpgs] = dnam_copy.loc[existing_cpgs].where(
            ~mask, means[existing_cpgs], axis=0
        )
    else:
        dnam_copy = dnam_copy.where(dnam_copy.notna(), means, axis=0)

    return dnam_copy


def hybrid_impute(dnam, cpg_source, required_cpgs, threshold=0.8):
    """
    Imputes missing values in a DNA methylation dataset based on a threshold. Sites with data below the threshold are replaced from an external source, while others are imputed using the average of existing values.

    Args:
        dnam (pd.DataFrame): DataFrame with samples as columns and CpG sites as rows.
        cpg_source (pd.Series): Series containing reference values for CpG sites.
        required_cpgs (list of str): List of CpG sites that need to be in the final dataset.
        threshold (float, optional): Threshold for determining imputation strategy. Default is 0.8.

    Returns:
        pd.DataFrame: DataFrame with missing values filled.

    Raises:
        ValueError: If certain required CpG sites are missing from both the dataset and the cpg_source.
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
    """
    Impute missing values in a DNA methylation dataset using a standard Biolearn dataset for reference.

    Args:
        dnam (pd.DataFrame): DataFrame with samples as columns and CpG sites as rows.

    Returns:
        pd.DataFrame: DataFrame with missing values filled.
    """
    biolearn_averages_file = get_data_file("biolearn_averages_450k.csv")
    df = pd.read_csv(biolearn_averages_file, index_col=0)
    biolearn_averages = df["average"]
    return impute_from_standard(dnam, biolearn_averages)
