import pandas as pd
import numpy as np
from biolearn.util import get_data_file
from biolearn.dunedin_pace import dunedin_pace_normalization


def anti_trafo(x, adult_age=20):
    y = np.where(
        x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age
    )
    return y


model_definitions = {
    "Horvathv1": {
        "year": 2013,
        "species": "Human",
        "tissue": "Multi-tissue",
        "source": "https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Horvath1.csv",
            "transform": lambda sum: anti_trafo(sum + 0.696),
        },
    },
    "Hannum": {
        "year": 2013,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.sciencedirect.com/science/article/pii/S1097276512008933",
        "output": "Age (Years)",
        "model": {"type": "LinearMethylationModel", "file": "Hannum.csv"},
    },
    "Lin": {
        "year": 2016,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.aging-us.com/article/100908/text",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Lin.csv",
            "transform": lambda sum: sum + 12.2169841,
        },
    },
    "PhenoAge": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.aging-us.com/article/101414/text",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "PhenoAge.csv",
            "transform": lambda sum: sum + 60.664,
        },
    },
    "Horvathv2": {
        "year": 2018,
        "species": "Human",
        "tissue": "Skin + blood",
        "source": "https://www.aging-us.com/article/101508/text",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Horvath2.csv",
            "transform": lambda sum: anti_trafo(sum - 0.447119319),
        },
    },
    "PEDBE": {
        "year": 2019,
        "species": "Human",
        "tissue": "Buccal",
        "source": "https://www.pnas.org/doi/10.1073/pnas.1820843116",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "PEDBE.csv",
            "transform": lambda sum: anti_trafo(sum - 2.1),
        },
    },
    "Zhang_10": {
        "year": 2019,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.nature.com/articles/ncomms14617",
        "output": "Mortality Risk",
        "model": {"type": "LinearMethylationModel", "file": "Zhang_10.csv"},
    },
    "DunedinPoAm38": {
        "year": 2020,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://elifesciences.org/articles/54870#s2",
        "output": "Aging Rate (Years/Year)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "DunedinPoAm38.csv",
            "transform": lambda sum: sum - 0.06929805,
        },
    },
    "DunedinPACE": {
        "year": 2022,
        "species": "Human",
        "tissue": "Unknown",
        "source": "https://www.proquest.com/docview/2634411178",
        "output": "Aging Rate (Years/Year)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "DunedinPACE.csv",
            "transform": lambda sum: sum - 1.949859,
            "preprocess": dunedin_pace_normalization,
        },
    },
    "AlcoholMcCartney": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6158884/",
        "output": "Alcohol Consumption",
        "model": {"type": "LinearMethylationModel", "file": "Alcohol.csv"},
    },
    "BMI_McCartney": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6158884/",
        "output": "BMI",
        "model": {"type": "LinearMethylationModel", "file": "BMI.csv"},
    },
    "DNAmTL": {
        "year": 2019,
        "species": "Human",
        "tissue": "Blood, Adipose",
        "source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6738410/",
        "output": "Telomere Length",
        "model": {
            "type": "LinearMethylationModel",
            "file": "DNAmTL.csv",
            "transform": lambda sum: sum - 7.924780053,
        },
    },
    "HRSInCHPhenoAge": {
        "year": "unknown",
        "species": "Human",
        "tissue": "unknown",
        "source": "unknown",
        "model": {
            "type": "LinearMethylationModel",
            "file": "HRSInCHPhenoAge.csv",
            "transform": lambda sum: sum + 52.8334080,
        },
    },
    "Knight": {
        "year": 2016,
        "species": "Human",
        "tissue": "Cord Blood",
        "source": "https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1068-z",
        "output": "Gestational Age",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Knight.csv",
            "transform": lambda sum: sum + 41.7,
        },
    },
    "LeeControl": {
        "year": 2019,
        "species": "Human",
        "tissue": "Placenta",
        "source": "https://www.aging-us.com/article/102049/text",
        "output": "Gestational Age",
        "model": {
            "type": "LinearMethylationModel",
            "file": "LeeControl.csv",
            "transform": lambda sum: sum + 13.06182,
        },
    },
    "LeeRefinedRobust": {
        "year": 2019,
        "species": "Human",
        "tissue": "Placenta",
        "source": "https://www.aging-us.com/article/102049/text",
        "output": "Gestational Age",
        "model": {
            "type": "LinearMethylationModel",
            "file": "LeeRefinedRobust.csv",
            "transform": lambda sum: sum + 30.74966,
        },
    },
    "LeeRobust": {
        "year": 2019,
        "species": "Human",
        "tissue": "Placenta",
        "source": "https://www.aging-us.com/article/102049/text",
        "output": "Gestational Age",
        "model": {
            "type": "LinearMethylationModel",
            "file": "LeeRobust.csv",
            "transform": lambda sum: sum + 24.99772,
        },
    },
    "SmokingMcCartney": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6158884/",
        "output": "Smoking Status",
        "model": {"type": "LinearMethylationModel", "file": "Smoking.csv"},
    },
}


class LinearMethylationModel:
    def __init__(
        self, coeffecient_file, transform, preprocess=None, **metadata
    ) -> None:
        self.transform = transform
        self.coefficients = pd.read_csv(
            get_data_file(coeffecient_file), index_col=0
        )
        self.preprocess = preprocess
        self.metadata = metadata

    @staticmethod
    def from_definition(clock_definition):
        def no_transform(_):
            return _

        model_def = clock_definition["model"]
        return LinearMethylationModel(
            model_def["file"],
            model_def.get("transform", no_transform),
            model_def.get("preprocess", no_transform),
            **{k: v for k, v in clock_definition.items() if k != "model"},
        )

    def predict(self, dnam_data):
        dnam_data = self.preprocess(dnam_data)

        # Join the coefficients and dnam_data on the index
        methylation_df = self.coefficients.join(dnam_data, how="inner")

        # Vectorized multiplication: multiply CoefficientTraining with all columns of dnam_data
        result = (
            methylation_df.iloc[:, 1:]
            .multiply(methylation_df["CoefficientTraining"], axis=0)
            .sum(axis=0)
        )
        return result.apply(self.transform)

    def methylation_sites(self):
        return list(self.coefficients.index)


class ImputationDecorator:
    def __init__(self, clock, imputation_method):
        self.clock = clock
        self.imputation_method = imputation_method

    def predict(self, dnam_data):
        # Impute missing values before prediction
        needed_cpgs = self.clock.methylation_sites()
        dnam_data_imputed = self.imputation_method(dnam_data, needed_cpgs)
        return self.clock.predict(dnam_data_imputed)

    # Forwarding other methods and attributes to the clock
    def __getattr__(self, name):
        return getattr(self.clock, name)


def single_sample_clock(clock_function, data):
    return clock_function(data).iloc[0, 0]
