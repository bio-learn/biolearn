import os
import cvxpy as cp
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


CLOCK_FOUNDATION_USAGE = "For cosmetics or life insurance applications, contact UCLA TDG regarding licensing status. For all other commercial usage `contact the Clock Foundation <https://clockfoundation.org/contact-us/>`_."

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
        },
        "usage": {
            "commercial": CLOCK_FOUNDATION_USAGE,
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
        "usage": {
            "commercial": CLOCK_FOUNDATION_USAGE,
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
        "usage": {
            "commercial": CLOCK_FOUNDATION_USAGE,
        },
    },
    "GrimAgeV2": {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9792204/",
        "output": "Mortality Adjusted Age (Years)",
        "model": {"type": "GrimageModel", "file": "GrimAgeV2.csv"},
        "usage": {
            "commercial": CLOCK_FOUNDATION_USAGE,
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
    "DNAmTL": {
        "year": 2019,
        "species": "Human",
        "tissue": "Blood, Adipose",
        "source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6738410/",
        "output": "Telomere Length",
        "model": {
            "type": "LinearMethylationModel",
            "file": "DNAmTL.csv",
        },
        "usage": {
            "commercial": CLOCK_FOUNDATION_USAGE,
        },
    },
    "HRSInCHPhenoAge": {
        "year": 2022,
        "species": "Human",
        "tissue": "Blood",
        "output": "Age (Years)",
        "source": "https://www.nature.com/articles/s43587-022-00248-2",
        "model": {
            "type": "LinearMethylationModel",
            "file": "HRSInCHPhenoAge.csv",
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
        },
    },
    "SexEstimation": {
        "year": 2021,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07675-2",
        "output": "Sex",
        "model": {
            "type": "SexEstimationModel",
            "file": "estimateSex.csv",
            "default_imputation": "averaging",
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
    "DeconvoluteBlood450K": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.biorxiv.org/content/10.1101/2023.12.02.569722v6",
        "output": "Cell Proportions",
        "model": {
            "type": "DeconvolutionModel",
            "file": "450K_reinius_12_reference.csv",
            "platform": "450K",
            "default_imputation": "none",
        },
        "usage": {
            "commercial": "Free to use",
            "non-commercial": "Free to use",
        },
    },
    "DeconvoluteBloodEPIC": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.biorxiv.org/content/10.1101/2023.12.02.569722v6",
        "output": "Cell Proportions",
        "model": {
            "type": "DeconvolutionModel",
            "file": "EPIC_salas_18_reference.csv",
            "platform": "EPIC",
            "default_imputation": "none",
        },
        "usage": {
            "commercial": "Free to use",
            "non-commercial": "Free to use",
        },
    },
    "TwelveCellDeconvoluteBloodEPIC": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.biorxiv.org/content/10.1101/2023.12.02.569722v6",
        "output": "Cell Proportions",
        "model": {
            "type": "DeconvolutionModel",
            "file": "twelve_cell_deconv.csv",
            "platform": "EPIC",
            "default_imputation": "none",
        },
        "usage": {
            "commercial": "Free to use",
            "non-commercial": "Free to use",
        },
    },
    "BMI_McCartney": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/s13059-018-1514-1",
        "output": "BMI",
        "model": {
            "type": "LinearMethylationModel",
            "file": "BMI_McCartney.csv",
            "transform": lambda x: 1 / (1 + np.exp(-x)),
        },
    },
    "EducationMcCartney": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/s13059-018-1514-1",
        "output": "Educational Attainment",
        "model": {
            "type": "LinearMethylationModel",
            "file": "EducationMcCartney.csv",
            "transform": lambda x: 1 / (1 + np.exp(-x)),
        },
    },
    "TotalCholesterolMcCartney": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/s13059-018-1514-1",
        "output": "Total Cholesterol",
        "model": {
            "type": "LinearMethylationModel",
            "file": "TotalCholesterolMcCartney.csv",
            "transform": lambda x: 1 / (1 + np.exp(-x)),
        },
    },
    "HDLCholesterolMcCartney": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/s13059-018-1514-1",
        "output": "HDL Cholesterol",
        "model": {
            "type": "LinearMethylationModel",
            "file": "HDLCholesterolMcCartney.csv",
            "transform": lambda x: 1 / (1 + np.exp(-x)),
        },
    },
    "LDLCholesterolMcCartney": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/s13059-018-1514-1",
        "output": "LDL with Remnant Cholesterol",
        "model": {
            "type": "LinearMethylationModel",
            "file": "LDLCholesterolMcCartney.csv",
            "transform": lambda x: 1 / (1 + np.exp(-x)),
        },
    },
    "BodyFatMcCartney": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/s13059-018-1514-1",
        "output": "Percentage Body Fat",
        "model": {
            "type": "LinearMethylationModel",
            "file": "BodyFatMcCartney.csv",
            "transform": lambda x: 1 / (1 + np.exp(-x)),
        },
    },
    "BMI_Reed": {
        "year": 2020,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/s13148-020-00841-5",
        "output": "BMI",
        "model": {
            "type": "LinearMethylationModel",
            "file": "BMI_Reed.csv",
        },
    },
    "ProstateCancerKirby": {
        "year": 2017,
        "species": "Human",
        "tissue": "Prostate",
        "source": "https://doi.org/10.1186/s12885-017-3252-2",
        "output": "Prostate Cancer Status",
        "model": {
            "type": "LinearMethylationModel",
            "file": "ProstateCancerKirby.csv",
        },
    },
    "HepatoXu": {
        "year": 2017,
        "species": "Human",
        "tissue": "Circulating DNA",
        "source": "https://doi.org/10.1038/nmat4997",
        "output": "Hepatocellular Carcinoma Status",
        "model": {
            "type": "LinearMethylationModel",
            "file": "HepatoXu.csv",
        },
    },
    "CVD_Westerman": {
        "year": 2020,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1161/JAHA.119.015299",
        "output": "Coronary Heart Disease Status",
        "model": {
            "type": "LinearMethylationModel",
            "file": "CVD_Westermann.csv",
            "transform": lambda x: 1 / (1 + np.exp(-x)),
        },
    },
    "AD_Bahado-Singh": {
        "year": 2021,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1371/journal.pone.0248375",
        "output": "Alzheimer's Disease Status",
        "model": {
            "type": "LinearMethylationModel",
            "file": "AD_Bahado-Singh.csv",
            "transform": lambda x: 1 / (1 + np.exp(-x - 0.072)),
        },
    },
    "DepressionBarbu": {
        "year": 2021,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1038/s41380-020-0808-3",
        "output": "Depression Risk",
        "model": {
            "type": "LinearMethylationModel",
            "file": "DepressionBarbu.csv",
        },
    },
    "TranscriptomicPredictionModel": {
        "year": 2015,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.nature.com/articles/ncomms9570",
        "output": "Age (Years)",
        "model": {
            "type": "LinearTranscriptomicModel",
            "file": "TranscriptomicPrediction.csv",
            "preprocess": lambda rna_data: preprocess_rna(
                map_ensembl_to_gene(rna_data)
            ),
            "transform": lambda sum: sum + 55.808884324,
        },
    },
}


def quantile_normalize(df):
    rank_mean = (
        df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
    )
    return df.rank(method="min").stack().astype(int).map(rank_mean).unstack()


def preprocess_rna(rna_matrix):
    normalized_data = quantile_normalize(rna_matrix)

    # Log2 transform the normalized data
    log2_data = np.log2(normalized_data + 1)  # Adding 1 to avoid log(0) error

    # Center the probe means to zero
    probe_centered_data = log2_data.sub(log2_data.mean(axis=1), axis="index")

    # Center the sample means to zero
    sample_centered_data = probe_centered_data.sub(
        probe_centered_data.mean(axis=0), axis="columns"
    )
    return sample_centered_data


def map_ensembl_to_gene(rna_matrix):
    mapping_file = get_data_file("reference/ensembl_to_gene.csv")
    ensembl_to_gene = pd.read_csv(mapping_file, index_col=0)
    id_to_gene_map = ensembl_to_gene["gene"].to_dict()
    new_rna_matrix = rna_matrix.rename(index=id_to_gene_map)
    return new_rna_matrix


class DeconvolutionModel:
    def __init__(self, reference_file, platform_input):

        # sets reference data
        self.reference = pd.read_csv(
            get_data_file(reference_file), index_col=0
        )

        # sets reference platform type (e.g. "450K")
        self.platform = platform_input

    @classmethod
    def from_definition(cls, clock_definition):

        model_def = clock_definition["model"]

        # calls constructor: "__init__()"
        return cls(model_def["file"], model_def["platform"])

    def predict(self, geo_data):

        # This function computes estimate of cell proportions using quadratic programming
        # inputs:
        # meth_vector: ndarray vector of bulk methylation (length must equal number of rows in deconv_reference)
        # deconv_reference: ndarray matrix with rows corresponding to CpG sites and columns to cell types
        # outputs:
        # cell_prop_estimate.value: ndarray vector of estimated cell type proportions
        def solve_qp(meth_vector, deconv_reference):

            # define cell proportion variable being solved for
            cell_prop_estimate = cp.Variable(deconv_reference.shape[1])

            # define quadratic objective function
            objective = cp.Minimize(
                cp.norm(deconv_reference @ cell_prop_estimate - meth_vector, 2)
            )

            # define constraints: proportion estimate should be within [0, 1] interval, sum to 1
            constraints = [
                0 <= cell_prop_estimate,
                cell_prop_estimate <= 1,
                cp.sum(cell_prop_estimate) == 1,
            ]

            # create optimization problem
            problem = cp.Problem(objective, constraints)

            # solve optimization problem
            problem.solve(solver=cp.ECOS)

            # return result
            return cell_prop_estimate.value

        # create copy of methylation matrix
        meth = geo_data.dnam.copy()

        # check if missing values are present in methylation matrix
        if any(meth.isna().values.flatten()):
            print(
                "WARNING: methylation data contains missing values. Rows with missing values will be removed"
            )

        # omit rows with missing values
        meth.dropna(inplace=True)

        # filter methylation array and reference to include only shared CpGs
        intersecting_cpgs = meth.index.intersection(self.reference.index)

        meth_filtered = meth.loc[intersecting_cpgs]

        reference_filtered = self.reference.loc[intersecting_cpgs]

        # confirm same number of rows in reference and methylation samples
        if meth_filtered.shape[0] != reference_filtered.shape[0]:
            diff = abs(meth_filtered.shape[0] - reference_filtered.shape[0])
            if meth_filtered.shape[0] > reference_filtered.shape[0]:
                bigger = "meth_filtered"
                smaller = "reference_filtered"
            else:
                bigger = "reference_filtered"
                smaller = "meth_filtered"

            msg = (
                f"Mismatch in row numbers: {bigger} has more rows than {smaller} "
                f"by {diff} rows."
            )

            if diff < 50:
                index_difference = set(meth_filtered.index) ^ set(
                    reference_filtered.index
                )
                msg += f"\nIndex differences (limited to 50): {list(index_difference)}"

            raise ValueError(msg)

        # save reference cell types
        cell_types = reference_filtered.columns

        # save methylation sample names
        sample_names = meth_filtered.columns

        # convert methylation samples and reference to ndarray
        meth_filtered_ndarray = meth_filtered.to_numpy()
        reference_filtered_ndarray = reference_filtered.to_numpy()

        # deconvolve each sample (each column of meth)
        cell_prop = np.apply_along_axis(
            solve_qp,
            axis=0,
            arr=meth_filtered_ndarray,
            deconv_reference=reference_filtered_ndarray,
        )

        # convert cell proportion ndarray to dataframe
        cell_prop_df = pd.DataFrame(cell_prop, columns=sample_names)
        cell_prop_df.index = cell_types

        return cell_prop_df

    # returns required methylation sites
    def methylation_sites(self):
        return list(self.reference.index)


class LinearModel:
    def __init__(
        self, coefficient_file_or_df, transform, preprocess=None, **metadata
    ) -> None:
        self.transform = transform
        self.preprocess = preprocess
        self.metadata = metadata

        if isinstance(coefficient_file_or_df, pd.DataFrame):
            self.coefficients = coefficient_file_or_df
        elif isinstance(coefficient_file_or_df, str):
            if os.path.isabs(coefficient_file_or_df):
                # Absolute path: load the file directly
                file_path = coefficient_file_or_df
            else:
                # Relative path: get the full path from the data folder
                file_path = get_data_file(coefficient_file_or_df)

            self.coefficients = pd.read_csv(file_path, index_col=0)
        else:
            raise ValueError(
                "coefficient_file_or_df must be either a DataFrame or a file path as a string."
            )

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
        matrix_data = self._get_data_matrix(geo_data)
        matrix_data = self.preprocess(matrix_data)
        matrix_data.loc["intercept"] = 1

        # Join the coefficients and dnam_data on the index
        model_df = self.coefficients.join(matrix_data, how="inner")

        # Vectorized multiplication: multiply CoefficientTraining with all columns of dnam_data
        result = (
            model_df.iloc[:, 1:]
            .multiply(model_df["CoefficientTraining"], axis=0)
            .sum(axis=0)
        )

        # Return as a DataFrame
        return result.apply(self.transform).to_frame(name="Predicted")

    def _get_data_matrix(self, geo_data):
        raise NotImplementedError()


class LinearMethylationModel(LinearModel):
    def _get_data_matrix(self, geo_data):
        return geo_data.dnam

    def methylation_sites(self):
        unique_vars = set(self.coefficients.index) - {"intercept"}
        return list(unique_vars)


class LinearTranscriptomicModel(LinearModel):
    def _get_data_matrix(self, geo_data):
        return geo_data.rna


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
