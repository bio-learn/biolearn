import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from biolearn.imputation import hybrid_impute
from biolearn.util import get_data_file


def dunedin_pace_normalization(dataframe):
    gold_standard_df = pd.read_csv(
        get_data_file("DunedinPACE_Gold_Means.csv"), index_col=0
    )
    gold_standard_means = dict(
        zip(gold_standard_df.index, gold_standard_df["mean"])
    )

    dataframe_transposed = dunedin_pace_preprocess_data(
        dataframe, gold_standard_df
    )

    # Keep only rows present in gold_standard_means
    dataframe_transposed = dataframe_transposed[
        dataframe_transposed.index.isin(gold_standard_means)
    ]

    # Ensure target is ordered according to dataframe_transposed's rows
    target_ordered = dataframe_transposed.index.map(
        gold_standard_means
    ).tolist()
    # Store the original index and columns
    original_index = dataframe_transposed.index
    original_columns = dataframe_transposed.columns

    dataframe_normalized = quantile_normalize_using_target(
        dataframe_transposed.values, target_ordered
    )
    dataframe_normalized = pd.DataFrame(
        dataframe_normalized,
        columns=original_columns,
        index=original_index,
    )

    return dataframe_normalized


def quantile_normalize_using_target(data, target_values):
    """
    Apply quantile normalization on data using target values.
    """
    sorted_target = np.sort(target_values)

    for _, column_data in enumerate(data.T):
        ranks = rankdata(column_data)

        rank_floor_values = np.floor(ranks).astype(int)
        has_decimal_above_0_4 = (ranks - rank_floor_values) > 0.4

        column_data[has_decimal_above_0_4] = 0.5 * (
            sorted_target[rank_floor_values[has_decimal_above_0_4] - 1]
            + sorted_target[rank_floor_values[has_decimal_above_0_4]]
        )
        column_data[~has_decimal_above_0_4] = sorted_target[
            rank_floor_values[~has_decimal_above_0_4] - 1
        ]

    return data


def dunedin_pace_preprocess_data(betas, means, proportionOfProbesRequired=0.8):
    if not (betas.map(np.isreal)).all().all():
        raise ValueError("betas matrix/data.frame is not numeric!")

    coefficients = pd.read_csv(get_data_file("DunedinPACE.csv"), index_col=0)
    model_probes = coefficients.index.tolist()
    gold_standard_probes = means.index.tolist()

    probeOverlap = len(set(betas.index) & set(model_probes)) / len(
        set(model_probes)
    )
    probeOverlap_background = len(
        set(betas.index) & set(gold_standard_probes)
    ) / len(set(gold_standard_probes))

    if (
        probeOverlap < proportionOfProbesRequired
        or probeOverlap_background < proportionOfProbesRequired
    ):

        warnings.warn(
            "Number of datapoints does not meet minimum required for this clock"
        )

    desired_labels = gold_standard_probes
    existing_labels = betas.index.intersection(desired_labels)
    betas_mat = betas.loc[existing_labels]
    return hybrid_impute(betas_mat, means["mean"], gold_standard_probes)
