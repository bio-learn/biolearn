import os
import pandas as pd


def horvath_transform(mult_sum):
    const = 0.695507258
    BA = (mult_sum + const) * 21 + 20
    return BA


def phenoage_transform(mult_sum):
    const = 60.664
    BA = mult_sum + const
    return BA


def no_transform(_):
    return _


def run_clock(dataframe, coeffecient_file, transform_function):
    script_dir = os.path.dirname(__file__)  # get the directory of the current script
    data_file_path = os.path.join(
        script_dir, "data", coeffecient_file
    )  # build the path to the data file

    coefficients = pd.read_csv(data_file_path, index_col=0)
    methylation_df = coefficients.merge(
        dataframe.transpose(), left_index=True, right_index=True
    )
    for c in methylation_df.columns[1:]:
        methylation_df[c] = methylation_df["CoefficientTraining"] * methylation_df[c]
    df_sum = methylation_df.drop("CoefficientTraining", axis=1).sum()
    return df_sum.apply(transform_function).to_frame(name="biological_age")


def horvath_clock(dataframe):
    """Runs the Horvath DNA methylation clock on each individual in the dataset to predict a biological age
    Parameters
    ----------
    dataframe : Pandas.Dataframe
        A pandas dataframe where each row represents an individual and each column represents a measurement about that individual.
        Needs to have DNA methylation measurements for the clock to work
    
    Returns
    -----------
    df: Pandas.Dataframe
        A pandas dataframe where each row represents an individual with a single column for predicted biological age
    """
    return run_clock(dataframe, "horvath.csv", horvath_transform)


def hannum_clock(dataframe):
    """Runs the Hannum DNA methylation clock on each individual in the dataset to predict a biological age
    Parameters
    ----------
    dataframe : Pandas.Dataframe
        A pandas dataframe where each row represents an individual and each column represents a measurement about that individual.
        Needs to have DNA methylation measurements for the clock to work
    
    Returns
    -----------
    df: Pandas.Dataframe
        A pandas dataframe where each row represents an individual with a single column for predicted biological age
    """
    return run_clock(dataframe, "hannum.csv", no_transform)


def phenoage_clock(dataframe):
    """Runs the PhenoAge DNA methylation clock on each individual in the dataset to predict a biological age
    Parameters
    ----------
    dataframe : Pandas.Dataframe
        A pandas dataframe where each row represents an individual and each column represents a measurement about that individual.
        Needs to have DNA methylation measurements for the clock to work
    
    Returns
    -----------
    df: Pandas.Dataframe
        A pandas dataframe where each row represents an individual with a single column for predicted biological age
    """
    return run_clock(dataframe, "phenoage.csv", phenoage_transform)


def single_sample_clock(clock_function, data):
    return clock_function(data).iloc[0, 0]
