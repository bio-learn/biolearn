import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

from biolearn.data_library import GeoData
from biolearn.dunedin_pace import dunedin_pace_normalization
from biolearn.util import get_data_file


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
    "YingCausAge": {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.biorxiv.org/content/10.1101/2022.10.07.511382v2",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "YingCausAge.csv",
            "transform": lambda sum: sum + 86.80816381,
        },
    },
    "YingDamAge": {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.biorxiv.org/content/10.1101/2022.10.07.511382v2",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "YingDamAge.csv",
            "transform": lambda sum: sum + 543.4315887,
        },
    },
    "YingAdaptAge": {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.biorxiv.org/content/10.1101/2022.10.07.511382v2",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "YingAdaptAge.csv",
            "transform": lambda sum: sum - 511.9742762,
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
        "tissue": "Blood",
        "source": "https://www.proquest.com/docview/2634411178",
        "output": "Aging Rate (Years/Year)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "DunedinPACE.csv",
            "transform": lambda sum: sum - 1.949859,
            "preprocess": dunedin_pace_normalization,
            "default_imputation": "none",
        },
    },
    "GrimAgeV1": {
        "year": 2019,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6366976/",
        "output": "Mortality Adjusted Age (Years)",
        "model": {"type": "GrimageModel", "file": "GrimAgeV1.csv"},
    },
    "GrimAgeV2": {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9792204/",
        "output": "Mortality Adjusted Age (Years)",
        "model": {"type": "GrimageModel", "file": "GrimAgeV2.csv"},
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
        "year": "2022",
        "species": "Human",
        "tissue": "Blood",
        "output": "Age (Years)",
        "source": "https://www.nature.com/articles/s43587-022-00248-2",
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
    "SexEstimation": {
        "year": 2021,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07675-2",
        "output": "Sex",
        "model": {"type": "SexEstimationModel", "file": "estimateSex.csv"},
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
    "DownSyndrome": {
        "year": 2021,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.nature.com/articles/s41467-021-21064-z",
        "output": "Down Syndrome Prediction",
        "model": {
            "type": "LinearMethylationModel",
            "file": "down_syndrome.csv",
            "default_imputation": "averaging",
        },
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

    @classmethod
    def from_definition(cls, clock_definition):
        def no_transform(_):
            return _

        model_def = clock_definition["model"]
        return cls(
            model_def["file"],
            model_def.get("transform", no_transform),
            model_def.get("preprocess", no_transform),
            **{k: v for k, v in clock_definition.items() if k != "model"},
        )

    def predict(self, geo_data):
        dnam_data = self.preprocess(geo_data.dnam)

        # Join the coefficients and dnam_data on the index
        methylation_df = self.coefficients.join(dnam_data, how="inner")

        # Vectorized multiplication: multiply CoefficientTraining with all columns of dnam_data
        result = (
            methylation_df.iloc[:, 1:]
            .multiply(methylation_df["CoefficientTraining"], axis=0)
            .sum(axis=0)
        )

        # Return as a DataFrame
        return result.apply(self.transform).to_frame(name="Predicted")

    def methylation_sites(self):
        return list(self.coefficients.index)


class GrimageModel:
    def __init__(self, coefficient_file, **metadata):
        self.coefficients = pd.read_csv(
            get_data_file(coefficient_file), index_col=0
        )
        self.metadata = metadata

    @classmethod
    def from_definition(cls, clock_definition):
        model_def = clock_definition["model"]
        return cls(
            model_def["file"],
            **{k: v for k, v in clock_definition.items() if k != "model"},
        )

    def predict(self, geo_data):
        if "sex" not in geo_data.metadata or "age" not in geo_data.metadata:
            raise ValueError("Metadata must contain 'sex' and 'age' columns")

        df = geo_data.dnam

        # Transposing metadata so that its structure aligns with dnam (columns as samples)
        transposed_metadata = geo_data.metadata.transpose()

        # Add metadata rows to dnam DataFrame
        df.loc["Age"] = transposed_metadata.loc["age"]
        df.loc["Female"] = transposed_metadata.loc["sex"].apply(
            lambda x: 1 if x == 1 else 0
        )
        df.loc["Intercept"] = 1

        grouped = self.coefficients.groupby("Y.pred")
        all_data = pd.DataFrame()

        for name, group in grouped:
            if name == "COX":
                cox_coefficients = group.set_index("var")["beta"]
                print(cox_coefficients)
            elif name == "transform":
                transform = group.set_index("var")["beta"]
                m_age = transform["m_age"]
                sd_age = transform["sd_age"]
                m_cox = transform["m_cox"]
                sd_cox = transform["sd_cox"]
            else:
                sub_clock_result = self.calculate_sub_clock(df, group)
                all_data[name] = sub_clock_result

        all_data["Age"] = geo_data.metadata["age"]
        all_data["Female"] = geo_data.metadata["sex"].apply(
            lambda x: 1 if x == 1 else 0
        )

        all_data["COX"] = all_data.mul(cox_coefficients).sum(axis=1)
        age_key = "DNAmGrimAge"
        accel_key = "AgeAccelGrim"
        # Calculate DNAmGrimAge
        Y = (all_data["COX"] - m_cox) / sd_cox
        all_data[age_key] = (Y * sd_age) + m_age

        # Calculate AgeAccelGrim
        lm = LinearRegression().fit(
            all_data[["Age"]].values, all_data[age_key].values
        )
        predictions = lm.predict(all_data[["Age"]].values)
        all_data[accel_key] = all_data[age_key] - predictions

        # Drop COX column after computations
        all_data.drop("COX", axis=1, inplace=True)

        return all_data

    def calculate_sub_clock(self, df, coefficients):
        # Filter coefficients for only those present in df
        relevant_coefficients = coefficients[
            coefficients["var"].isin(df.index)
        ]

        # Create a Series from the relevant coefficients, indexed by 'var'
        coefficients_series = relevant_coefficients.set_index("var")["beta"]

        # Align coefficients with df's rows and multiply, then sum across CpG sites for each sample
        result = (
            df.loc[coefficients_series.index]
            .multiply(coefficients_series, axis=0)
            .sum()
        )

        return result

    def rename_columns(self, data, old_names, new_names):
        for old_name, new_name in zip(old_names, new_names):
            data.rename(columns={old_name: new_name}, inplace=True)

    def methylation_sites(self):
        filtered_df = self.coefficients[
            ~self.coefficients.index.isin(["COX", "transform"])
        ]
        unique_vars = set(filtered_df["var"]) - {"Intercept", "Age", "Female"}
        return list(unique_vars)


class SexEstimationModel:
    def __init__(self, coeffecient_file, **metadata):
        self.coefficients = pd.read_csv(
            get_data_file(coeffecient_file), index_col=0, low_memory=False
        )
        self.metadata = metadata

    @classmethod
    def from_definition(cls, clock_definition):
        # Implementation for creating an instance from a definition
        # Adjust this as needed for your specific definition format
        return cls(
            clock_definition["model"]["file"],
            **{k: v for k, v in clock_definition.items() if k != "model"},
        )

    def predict(self, geo_data):
        dnam_data = geo_data.dnam
        # Find the common probes between the reference and the input data
        common_probes = dnam_data.index.intersection(self.coefficients.index)
        reference = self.coefficients.loc[common_probes]
        dnam_data = dnam_data.loc[common_probes]

        # Identify autosomes and calculate mean and standard deviation
        autosomes = reference.loc[~reference["CHR"].isin(["X", "Y"])].index
        d_mean = dnam_data.loc[autosomes].mean(axis=0, skipna=True)
        d_std = dnam_data.loc[autosomes].std(axis=0, skipna=True)

        # Normalize the data using Z-score normalization
        z_data = (
            dnam_data.subtract(d_mean, axis=1).div(d_std, axis=1).fillna(0)
        )

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

    def methylation_sites(self):
        return list(self.coefficients.index)


class ImputationDecorator:
    def __init__(self, clock, imputation_method):
        self.clock = clock
        self.imputation_method = imputation_method

    def predict(self, geo_data):
        # Impute missing values before prediction
        needed_cpgs = self.clock.methylation_sites()
        dnam_data_imputed = self.imputation_method(geo_data.dnam, needed_cpgs)

        return self.clock.predict(
            GeoData(geo_data.metadata, dnam_data_imputed)
        )

    # Forwarding other methods and attributes to the clock
    def __getattr__(self, name):
        return getattr(self.clock, name)


def single_sample_clock(clock_function, data):
    return clock_function(data).iloc[0, 0]
