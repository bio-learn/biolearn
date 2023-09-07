import os
from scipy.stats import rankdata
import pandas as pd
import numpy as np
from biolearn.util import get_data_file


def dunedin_pace_normalization(dataframe):
    gold_standard_df = pd.read_csv(get_data_file("DunedinPACE_Gold_Means.csv"))
    gold_standard_means = dict(
        zip(gold_standard_df["CpGmarker"], gold_standard_df["mean"])
    )

    original_index = dataframe.index
    dataframe_transposed = dunedin_pace_preprocess_data(
        dataframe.transpose(), gold_standard_df
    )

    # Keep only rows present in gold_standard_means
    dataframe_transposed = dataframe_transposed[
        dataframe_transposed.index.isin(gold_standard_means)
    ]

    # Ensure target is ordered according to dataframe_transposed's rows
    target_ordered = dataframe_transposed.index.map(
        gold_standard_means
    ).tolist()

    dataframe_normalized = quantile_normalize_using_target(
        dataframe_transposed.values, target_ordered
    )
    dataframe_normalized = pd.DataFrame(
        dataframe_normalized,
        columns=original_index,
        index=dataframe_transposed.index,
    )

    return dataframe_normalized.T


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

    gold_standard_means_lookup = dict(zip(means["CpGmarker"], means["mean"]))
    coefficients = pd.read_csv(get_data_file("DunedinPACE.csv"), index_col=0)
    model_probes = coefficients.index.tolist()
    gold_standard_probes = means["CpGmarker"].tolist()

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
        raise (
            Exception(
                "Number of datapoints does not meet minimum required for this clock"
            )
        )

    betas_mat = betas.loc[gold_standard_probes].copy()

    probesNotInMatrix = [
        probe for probe in gold_standard_probes if probe not in betas_mat.index
    ]

    for probe in probesNotInMatrix:
        betas_mat.loc[probe] = gold_standard_means_lookup[probe]

    samplesToRemove = betas_mat.columns[
        (betas_mat.notna().sum(axis=0) / len(betas_mat)).lt(
            proportionOfProbesRequired
        )
    ]

    betas_mat.drop(columns=samplesToRemove, inplace=True)

    pctValuesPresent = betas_mat.notna().mean(axis=1)

    probesToAdjust = pctValuesPresent.index[
        (pctValuesPresent < 1)
        & (pctValuesPresent >= proportionOfProbesRequired)
    ]
    for probe in probesToAdjust:
        betas_mat.loc[probe].fillna(betas_mat.loc[probe].mean(), inplace=True)

    probesToReplaceWithMean = pctValuesPresent.index[
        pctValuesPresent < proportionOfProbesRequired
    ]
    for probe in probesToReplaceWithMean:
        betas_mat.loc[probe] = gold_standard_means_lookup[probe]

    return betas_mat
