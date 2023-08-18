import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import quantile_transform

def horvath_transform(mult_sum):
    const = 0.695507258
    BA = (mult_sum + const) * 21 + 20
    return BA

def anti_trafo(x, adult_age=20):
    y = np.where(
        x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age
    )
    return y

def no_transform(_):
    return _

def get_data_file(filename):
    script_dir = os.path.dirname(
        __file__
    )  # get the directory of the current script
    return os.path.join(
        script_dir, "data", filename
    )  # build the path to the data file

def run_clock(dataframe, coeffecient_file, transform_function):
    coefficients = pd.read_csv(get_data_file(coeffecient_file), index_col=0)
    methylation_df = coefficients.merge(
        dataframe.transpose(), left_index=True, right_index=True
    )
    for c in methylation_df.columns[1:]:
        methylation_df[c] = (
            methylation_df["CoefficientTraining"] * methylation_df[c]
        )
    df_sum = methylation_df.drop("CoefficientTraining", axis=1).sum()
    return df_sum.apply(transform_function).to_frame(name="biological_age")


def horvathv1(dataframe):
    """Runs the Horvath DNA methylation clock on each individual in the dataset to predict a biological age

    Parameters
    ----------
    dataframe : Pandas.Dataframe
        A pandas dataframe where each row represents an individual and each column represents a measurement about that individual.
        Needs to have DNA methylation measurements for the clock to work

    Returns
    -------
    df: Pandas.Dataframe
        A pandas dataframe where each row represents an individual with a single column for predicted biological age
    """
    transform = lambda sum: anti_trafo(sum + 0.696)
    return run_clock(dataframe, "Horvath1.csv", transform)

def horvathv2(dataframe):
    transform = lambda sum: anti_trafo(sum - 0.447119319)
    return run_clock(dataframe, "Horvath2.csv", transform)

def hannum(dataframe):
    """Runs the Hannum DNA methylation clock on each individual in the dataset to predict a biological age

    Parameters
    ----------
    dataframe : Pandas.Dataframe
        A pandas dataframe where each row represents an individual and each column represents a measurement about that individual.
        Needs to have DNA methylation measurements for the clock to work

    Returns
    -------
    df: Pandas.Dataframe
        A pandas dataframe where each row represents an individual with a single column for predicted biological age
    """
    return run_clock(dataframe, "Hannum.csv", no_transform)


def phenoage(dataframe):
    """Runs the PhenoAge DNA methylation clock on each individual in the dataset to predict a biological age

    Parameters
    ----------
    dataframe : Pandas.Dataframe
        A pandas dataframe where each row represents an individual and each column represents a measurement about that individual.
        Needs to have DNA methylation measurements for the clock to work

    Returns
    -------
    df: Pandas.Dataframe
        A pandas dataframe where each row represents an individual with a single column for predicted biological age
    """
    transform = lambda sum: sum + 60.664
    return run_clock(dataframe, "PhenoAge.csv", transform)

# Results missing from expected file
# def bohlin(dataframe):
#     transform = lambda sum: sum + 277.2421
#     return run_clock(dataframe, "Bohlin.csv", transform)

def alcohol_mccartney(dataframe):
    return run_clock(dataframe, "Alcohol.csv", no_transform)

def bmi_mccartney(dataframe):
    return run_clock(dataframe, "BMI.csv", no_transform)

def dnam_tl(dataframe):
    transform = lambda sum: sum - 7.924780053
    return run_clock(dataframe, "DNAmTL.csv", transform)

import pandas as pd

def dunedin_pace_normalization(dataframe):
    # Read the gold standard means into a dictionary
    gold_standard_df = pd.read_csv(get_data_file("DunedinPACE_Means.csv"))
    gold_standard_means = dict(zip(gold_standard_df['CpGmarker'], gold_standard_df['mean']))
    
    # Identify cpgs in dataframe but not in gold standard, and drop them
    missing_cpgs = [cpg for cpg in dataframe.columns if cpg not in gold_standard_means]
    dataframe = dataframe.drop(columns=missing_cpgs)
    
    # Quantile Normalize the dataframe
    dataframe_normalized = quantile_transform(dataframe, axis=0, n_quantiles=dataframe.shape[0], copy=True, output_distribution='normal')
    dataframe_normalized = pd.DataFrame(dataframe_normalized, index=dataframe.index, columns=dataframe.columns)
    
    # Adjust the data to have the mean of the gold standard
    for cpg, mean_value in gold_standard_means.items():
        if cpg in dataframe_normalized.columns:
            dataframe_normalized[cpg] += mean_value - dataframe_normalized[cpg].mean()
    
    return dataframe_normalized

def dunedin_pace(dataframe):
    normalized_data = dunedin_pace_normalization(dataframe)
    transform = lambda sum: sum - 1.949859
    return run_clock(normalized_data, "DunedinPACE.csv", transform)

def dunedin_poam38(dataframe):
    transform = lambda sum: sum - 0.06929805
    return run_clock(dataframe, "DunedinPoAm38.csv", transform)

# Results missing from expected file
# def dnam_clock_cortical(dataframe):
#     transform = lambda sum: sum + 0.577682570446177
#     return run_clock(dataframe, "DNAmClockCortical.csv", transform)

def hrs_in_ch_phenoage(dataframe):
    transform = lambda sum: sum + 52.8334080
    return run_clock(dataframe, "HRSInCHPhenoAge.csv", transform)

def knight(dataframe):
    transform = lambda sum: sum + 41.7
    return run_clock(dataframe, "Knight.csv", transform)

def lee_control(dataframe):
    transform = lambda sum: sum + 13.06182
    return run_clock(dataframe, "LeeControl.csv", transform)

def lee_refined_robust(dataframe):
    transform = lambda sum: sum + 30.74966
    return run_clock(dataframe, "LeeRefinedRobust.csv", transform)

def lee_robust(dataframe):
    transform = lambda sum: sum + 24.99772
    return run_clock(dataframe, "LeeRobust.csv", transform)

def lin(dataframe):
    transform = lambda sum: sum + 12.2169841
    return run_clock(dataframe, "Lin.csv", transform)

# Test results do not match expected
# def mayne(dataframe):
#     transform = lambda sum: sum + 24.99026
#     return run_clock(dataframe, "Mayne.csv", transform)

# Coeffecient file is broken
# def mi_age(dataframe):
#     return run_clock(dataframe, "MiAge.csv", no_transform)

def pedbe(dataframe):
    transform = lambda sum: anti_trafo(sum - 2.1)
    return run_clock(dataframe, "PEDBE.csv", transform)

def phenoage(dataframe):
    """Runs the PhenoAge DNA methylation clock on each individual in the dataset to predict a biological age

    Parameters
    ----------
    dataframe : Pandas.Dataframe
        A pandas dataframe where each row represents an individual and each column represents a measurement about that individual.
        Needs to have DNA methylation measurements for the clock to work

    Returns
    -------
    df: Pandas.Dataframe
        A pandas dataframe where each row represents an individual with a single column for predicted biological age
    """
    transform = lambda sum: sum + 60.664
    return run_clock(dataframe, "PhenoAge.csv", transform)

def smoking_mccartney(dataframe):
    return run_clock(dataframe, "Smoking.csv", no_transform)

# Test results do not match expected
# def vidal_bralo(dataframe):
#     transform = lambda sum: sum + 84.7
#     return run_clock(dataframe, "Smoking.csv", transform)

def zhang_10(dataframe):
    return run_clock(dataframe, "Zhang_10.csv", no_transform)

# Test results do not match expected
# def zhang_2019(dataframe):
#     transform = lambda sum: sum + 65.8
#     return run_clock(dataframe, "Zhang2019.csv", transform)


def single_sample_clock(clock_function, data):
    return clock_function(data).iloc[0, 0]
