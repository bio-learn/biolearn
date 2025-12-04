import pandas as pd
import torch
from biolearn.util import get_data_file


def impute_from_standard(dnam, cpg_source, cpgs_to_impute=None):
    """
    Impute all missing values in a DNA methylation dataset using the reference values from an external dataset.

    Args:
        dnam (pd.DataFrame): DataFrame with samples as columns and CpG sites as rows.
        cpg_source (pd.Series): Series containing reference values for CpG sites.
        cpgs_to_impute (list of str, optional): List of CpG sites to impute. Missing cpgs will only be imputed if in this list.

    Returns:
        pd.DataFrame: DataFrame with missing values filled.
    """
    if cpgs_to_impute:
        dnam = dnam.reindex(dnam.index.union(cpgs_to_impute))
        impute_rows = dnam.loc[cpgs_to_impute]
        impute_rows = impute_rows.apply(lambda col: col.fillna(cpg_source))
        df_filled = dnam.combine_first(impute_rows)
    else:
        df_filled = dnam.apply(lambda col: col.fillna(cpg_source))
    return df_filled


def impute_from_average(dnam, cpgs_to_impute=None):
    """
    Impute all missing values in a DNA methylation dataset using the average from the dataset itself.

    Args:
        dnam (pd.DataFrame): DataFrame with samples as columns and CpG sites as rows.
        cpgs_to_impute (list of str, optional): List of CpG sites to impute. Missing cpgs will not be imputed.

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
        required_cpgs (list of str): List of CpG sites to impute. Missing cpgs will only be imputed if in this list.
        threshold (float, optional): Threshold for determining imputation strategy. Default is 0.8.

    Returns:
        pd.DataFrame: DataFrame with missing values filled.

    Raises:
        ValueError: If certain required CpG sites are missing from both the dataset and the cpg_source.
    """
    # Filter out rows below the threshold
    mask_above_threshold = dnam.notna().mean(axis=1) >= threshold
    dnam_filtered = dnam[mask_above_threshold]

    # Impute remaining rows using impute_from_average
    df_filled = impute_from_average(dnam_filtered)

    missing_cpgs_from_dataset = set(required_cpgs) - set(df_filled.index)
    missing_cpgs_from_source = [
        cpg for cpg in missing_cpgs_from_dataset if cpg not in cpg_source
    ]

    if missing_cpgs_from_source:
        raise ValueError(
            f"Tried to fill the following cpgs but they were missing from cpg_source: {missing_cpgs_from_source}"
        )

    # Create a DataFrame for missing CpGs with the same columns as dnam
    missing_cpgs_to_fill = {
        cpg: [cpg_source[cpg]] * len(dnam.columns)
        for cpg in missing_cpgs_from_dataset
    }
    missing_cpgs_df = pd.DataFrame.from_dict(
        missing_cpgs_to_fill, orient="index", columns=dnam.columns
    )

    # Exclude empty or all-NA columns before concatenation
    missing_cpgs_df = missing_cpgs_df.dropna(how="all", axis=1)

    # Concatenate the filled DataFrame with the missing CpGs DataFrame
    df_filled = pd.concat([df_filled, missing_cpgs_df])

    return df_filled.sort_index()


def impute_with_center(dnam, center, required_cpgs):
    """
    Imputes missing CpG sites using the center reference values.

    Args:
        dnam (DataFrame): The DNA methylation data with CpG sites as rows and samples as columns.
        center (Tensor or Series): A PyTorch tensor or pandas Series containing the center reference values for CpG sites.
        required_cpgs (list): A list of CpG sites required for the model.

    Returns:
        DataFrame: The DNA methylation data with missing CpG sites imputed.
    """
    # Convert center to a pandas Series if it's a tensor
    if isinstance(center, torch.Tensor):
        center = pd.Series(center.numpy(), index=required_cpgs)

    # Debug: Check the input data
    print("Input data shape:", dnam.shape)
    print("Input CpG sites (first 10):", dnam.index[:10])

    # Identify missing CpG sites in the dataset
    missing_cpgs_from_dataset = [
        cpg for cpg in required_cpgs if cpg not in dnam.index
    ]

    # Debug: Check missing CpGs from the dataset
    print(
        "Missing CpGs from dataset (first 10):", missing_cpgs_from_dataset[:10]
    )
    print(f"Total missing CpGs from dataset: {len(missing_cpgs_from_dataset)}")

    # Ensure all required CpG sites are in the center reference
    missing_cpgs_from_center = [
        cpg for cpg in missing_cpgs_from_dataset if cpg not in center.index
    ]

    if missing_cpgs_from_center:
        raise ValueError(
            f"Tried to fill the following CpGs but they were missing from the center reference: {missing_cpgs_from_center}"
        )

    # Handle the case where all CpGs are missing
    if len(missing_cpgs_from_dataset) == len(required_cpgs):
        print(
            "All required CpGs are missing from the dataset. Imputing all CpGs."
        )
        # Create a DataFrame for all required CpGs using the center values
        dnam_filled = pd.DataFrame(
            {cpg: [center[cpg]] * len(dnam.columns) for cpg in required_cpgs},
            index=required_cpgs,
            columns=dnam.columns,
        )
        return dnam_filled

    # Create a DataFrame for missing CpGs with the same columns as dnam
    missing_cpgs_to_fill = {
        cpg: [center[cpg]] * len(dnam.columns)
        for cpg in missing_cpgs_from_dataset
    }
    missing_cpgs_df = pd.DataFrame(missing_cpgs_to_fill, index=dnam.columns).T

    # Set the correct index for missing CpGs
    missing_cpgs_df.index = missing_cpgs_from_dataset

    # Append the missing CpGs to the original dnam DataFrame
    dnam_filled = pd.concat([dnam, missing_cpgs_df])

    # Reorder rows to match the required CpG order
    dnam_filled = dnam_filled.loc[required_cpgs]

    return dnam_filled
