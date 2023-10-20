import pandas as pd
import numpy as np
from biolearn.util import get_data_file
from biolearn.dunedin_pace import dunedin_pace_normalization


def anti_trafo(x, adult_age=20):
    y = np.where(
        x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age
    )
    return y


clock_definitions = {
    "Horvathv1": {
        "year": 2013,
        "species": "Human",
        "tissue": "Multi-tissue",
        "source": "https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115",
        "model": {
            "type": "LinearMethylationClock",
            "file": "Horvath1.csv",
            "transform": lambda sum: anti_trafo(sum + 0.696),
        },
    },
    "Hannum": {
        "year": 2013,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.sciencedirect.com/science/article/pii/S1097276512008933",
        "model": {"type": "LinearMethylationClock", "file": "Hannum.csv"},
    },
    "Lin": {
        "year": 2016,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.aging-us.com/article/100908/text",
        "model": {
            "type": "LinearMethylationClock",
            "file": "Lin.csv",
            "transform": lambda sum: sum + 12.2169841,
        },
    },
    "PhenoAge": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.aging-us.com/article/101414/text",
        "model": {
            "type": "LinearMethylationClock",
            "file": "PhenoAge.csv",
            "transform": lambda sum: sum + 60.664,
        },
    },
    "Horvathv2": {
        "year": 2018,
        "species": "Human",
        "tissue": "Skin + blood",
        "source": "https://www.aging-us.com/article/101508/text",
        "model": {
            "type": "LinearMethylationClock",
            "file": "Horvath2.csv",
            "transform": lambda sum: anti_trafo(sum - 0.447119319),
        },
    },
    "PEDBE": {
        "year": 2019,
        "species": "Human",
        "tissue": "Buccal",
        "source": "https://www.pnas.org/doi/10.1073/pnas.1820843116",
        "model": {
            "type": "LinearMethylationClock",
            "file": "PEDBE.csv",
            "transform": lambda sum: anti_trafo(sum - 2.1),
        },
    },
    "Zhang_10": {
        "year": 2019,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0667-1",
        "model": {"type": "LinearMethylationClock", "file": "Zhang_10.csv"},
    },
    "DunedinPoAm38": {
        "year": 2020,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://elifesciences.org/articles/54870#s2",
        "model": {
            "type": "LinearMethylationClock",
            "file": "DunedinPoAm38.csv",
            "transform": lambda sum: sum - 0.06929805,
        },
    },
    "DunedinPACE": {
        "year": 2022,
        "species": "Human",
        "tissue": "Unknown",
        "source": "https://www.proquest.com/docview/2634411178",
        "model": {
            "type": "LinearMethylationClock",
            "file": "DunedinPACE.csv",
            "transform": lambda sum: sum - 1.949859,
            "preprocess": dunedin_pace_normalization,
        },
    },
    "AlcoholMcCartney": {
        "year": "unknown",
        "species": "Human",
        "tissue": "unknown",
        "source": "unknown",
        "model": {"type": "LinearMethylationClock", "file": "Alcohol.csv"},
    },
    "BMI_McCartney": {
        "year": "unknown",
        "species": "Human",
        "tissue": "unknown",
        "source": "unknown",
        "model": {"type": "LinearMethylationClock", "file": "BMI.csv"},
    },
    "DNAmTL": {
        "year": "unknown",
        "species": "Human",
        "tissue": "unknown",
        "source": "unknown",
        "model": {
            "type": "LinearMethylationClock",
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
            "type": "LinearMethylationClock",
            "file": "HRSInCHPhenoAge.csv",
            "transform": lambda sum: sum + 52.8334080,
        },
    },
    "Knight": {
        "year": "unknown",
        "species": "Human",
        "tissue": "unknown",
        "source": "unknown",
        "model": {
            "type": "LinearMethylationClock",
            "file": "Knight.csv",
            "transform": lambda sum: sum + 41.7,
        },
    },
    "LeeControl": {
        "year": "unknown",
        "species": "Human",
        "tissue": "unknown",
        "source": "unknown",
        "model": {
            "type": "LinearMethylationClock",
            "file": "LeeControl.csv",
            "transform": lambda sum: sum + 13.06182,
        },
    },
    "LeeRefinedRobust": {
        "year": "unknown",
        "species": "Human",
        "tissue": "unknown",
        "source": "unknown",
        "model": {
            "type": "LinearMethylationClock",
            "file": "LeeRefinedRobust.csv",
            "transform": lambda sum: sum + 30.74966,
        },
    },
    "LeeRobust": {
        "year": "unknown",
        "species": "Human",
        "tissue": "unknown",
        "source": "unknown",
        "model": {
            "type": "LinearMethylationClock",
            "file": "LeeRobust.csv",
            "transform": lambda sum: sum + 24.99772,
        },
    },
    "SmokingMcCartney": {
        "year": "unknown",
        "species": "Human",
        "tissue": "unknown",
        "source": "unknown",
        "model": {"type": "LinearMethylationClock", "file": "Smoking.csv"},
    },
}


class LinearMethylationClock:
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
        return LinearMethylationClock(
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


def single_sample_clock(clock_function, data):
    return clock_function(data).iloc[0, 0]
