import pandas as pd
import numpy as np
from biolearn.util import get_data_file
from biolearn.dunedin_pace import dunedin_pace_normalization


def anti_trafo(x, adult_age=20):
    y = np.where(
        x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age
    )
    return y

def no_transform(_):
    return _

clock_definitions = {
    "Horvathv1": {"file": "Horvath1.csv", "transform": lambda sum: anti_trafo(sum + 0.696)},
    "Horvathv2": {"file": "Horvath2.csv", "transform": lambda sum: anti_trafo(sum - 0.447119319)},
    "Hannum": {"file": "Hannum.csv", "transform": no_transform},
    "PhenoAge": {"file": "PhenoAge.csv", "transform": lambda sum: sum + 60.664},
    "AlcoholMcCartney": {"file": "Alcohol.csv", "transform": no_transform},
    "BMI_McCartney": {"file": "BMI.csv", "transform": no_transform},
    "DNAmTL": {"file": "DNAmTL.csv", "transform": lambda sum: sum - 7.924780053},
    "DunedinPACE": {"file": "DunedinPACE.csv", "transform": lambda sum: sum - 1.949859, "preprocess": dunedin_pace_normalization},
    "DunedinPoAm38": {"file": "DunedinPoAm38.csv", "transform": lambda sum: sum - 0.06929805},
    "HRSInCHPhenoAge": {"file": "HRSInCHPhenoAge.csv", "transform": lambda sum: sum + 52.8334080},
    "Knight": {"file": "Knight.csv", "transform": lambda sum: sum + 41.7},
    "LeeControl": {"file": "LeeControl.csv", "transform": lambda sum: sum + 13.06182},
    "LeeRefinedRobust": {"file": "LeeRefinedRobust.csv", "transform": lambda sum: sum + 30.74966},
    "LeeRobust": {"file": "LeeRobust.csv", "transform": lambda sum: sum + 24.99772},
    "Lin": {"file": "Lin.csv", "transform": lambda sum: sum + 12.2169841},
    "PEDBE": {"file": "PEDBE.csv", "transform": lambda sum: anti_trafo(sum - 2.1)},
    "SmokingMcCartney": {"file": "Smoking.csv", "transform": no_transform},
    "Zhang_10": {"file": "Zhang_10.csv", "transform": no_transform},
}

class LinearMethylationClock:
    def __init__(self, coeffecient_file, transform, preprocess=None) -> None:
        self.transform = transform
        self.coefficients = pd.read_csv(coeffecient_file, index_col=0)
        self.preprocess = preprocess

    def predict(self, dnam_data):
        if self.preprocess:
            dnam_data = self.preprocess(dnam_data)

        # Join the coefficients and dnam_data on the index
        methylation_df = self.coefficients.join(dnam_data, how='inner')

        # Vectorized multiplication: multiply CoefficientTraining with all columns of dnam_data
        result = methylation_df.iloc[:, 1:].multiply(methylation_df["CoefficientTraining"], axis=0).sum(axis=0)
        return result.apply(self.transform)

    def methylation_sites(self):
        return list(self.coefficients.index)


def single_sample_clock(clock_function, data):
    return clock_function(data).iloc[0, 0]
