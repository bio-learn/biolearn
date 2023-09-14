from biolearn.util import get_data_file
import pandas as pd
import numpy as np


def estimate_sex(dnam_data):
    # Load the reference data
    reference = pd.read_csv(
        get_data_file("estimateSex.csv"), index_col=0, low_memory=False
    )

    # Find the common probes between the reference and the input data
    common_probes = dnam_data.index.intersection(reference.index)
    reference = reference.loc[common_probes]
    dnam_data = dnam_data.loc[common_probes]

    # Identify autosomes and calculate mean and standard deviation
    autosomes = reference.loc[~reference["CHR"].isin(["X", "Y"])].index
    d_mean = dnam_data.loc[autosomes].mean(axis=0, skipna=True)
    d_std = dnam_data.loc[autosomes].std(axis=0, skipna=True)

    # Normalize the data using Z-score normalization
    z_data = dnam_data.subtract(d_mean, axis=1).div(d_std, axis=1).fillna(0)

    # Perform the sex prediction for chromosomes X and Y separately
    pred_xy = {}
    for chr in ["X", "Y"]:
        chr_ref = reference.loc[reference["pca"] == chr]
        pred = (z_data.loc[chr_ref.index].T - chr_ref["mean"].values).dot(
            chr_ref["coeff"].values
        )
        pred_xy[chr] = pred

    # Create a DataFrame from the prediction results
    pred_df = pd.DataFrame(pred_xy)

    # Map the prediction results to sex categories
    pred_df["predicted_sex"] = "Female"
    pred_df.loc[
        (pred_df["X"] < 0) & (pred_df["Y"] > 0), "predicted_sex"
    ] = "Male"
    pred_df.loc[
        (pred_df["X"] > 0) & (pred_df["Y"] > 0), "predicted_sex"
    ] = "47,XXY"
    pred_df.loc[
        (pred_df["X"] < 0) & (pred_df["Y"] < 0), "predicted_sex"
    ] = "45,XO"

    return pred_df
