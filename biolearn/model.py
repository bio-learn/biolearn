import os
import cvxpy as cp
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.optimize import minimize_scalar
from sklearn.linear_model import LinearRegression
import requests
import io
import warnings
from typing import Optional, Dict, Any, Union


from biolearn.data_library import GeoData
from biolearn.dunedin_pace import dunedin_pace_normalization
from biolearn.util import get_data_file


def anti_trafo(x, adult_age=20):
    y = np.where(
        x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age
    )
    return y


def preprocess_pasta(df):
    df = df.loc[~df.index.duplicated(keep="first")]  # keep first duplicate
    genes = pd.read_csv("biolearn/data/Pasta.csv", index_col=0).index
    out = df.reindex(genes)  # select and order genes
    med = np.nanmedian(out.to_numpy())
    out = out.fillna(med)  # fill NA by overall median
    out = out.rank(
        axis=0, method="average", na_option="keep", ascending=True
    )  # rank normalization
    return out


def olink_standardization_preprocess(reference_sd_file):
    """
    Create a preprocessing function that standardizes Olink protein data
    to match reference standard deviations.

    Parameters:
    -----------
    reference_sd_file : str
        Path to CSV file containing reference standard deviations.
        File should have 'protein' and 'sd' columns.

    Returns:
    --------
    preprocess : function
        A preprocessing function that takes a DataFrame and returns
        a standardized DataFrame.
    """

    def preprocess(df):
        reference_path = get_data_file(reference_sd_file)
        reference_df = pd.read_csv(reference_path)
        reference_sds = reference_df.set_index("protein")["sd"]

        # Make a copy to avoid modifying the original
        df_standardized = df.copy()

        for col in df.columns:
            if col in reference_sds.index:
                current_sd = df[col].std(ddof=1, skipna=True)

                if current_sd > 0:  # Avoid division by zero
                    df_standardized[col] = (
                        df[col] / current_sd
                    ) * reference_sds[col]
            else:
                df_standardized[col] = np.nan

        return df_standardized

    return preprocess


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
    "GPAge10": {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00211-4",
        "output": "Age (Years)",
        "model": {
            "type": "GPAgeModel",
            "file": "GP-age_model_10_cpgs.json.zip",
            "sites": [
                "cg16867657",
                "cg06639320",
                "cg04875128",
                "cg19283806",
                "cg07553761",
                "cg08128734",
                "cg12934382",
                "cg00573770",
                "cg23479922",
                "cg10501210",
            ],
            "default_imputation": "none",
        },
    },
    "GPAge30": {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00211-4",
        "output": "Age (Years)",
        "model": {
            "type": "GPAgeModel",
            "file": "GP-age_model_30_cpgs.json.zip",
            "sites": [
                "cg16867657",
                "cg22454769",
                "cg06639320",
                "cg04875128",
                "cg19283806",
                "cg24724428",
                "cg07553761",
                "cg24079702",
                "cg08128734",
                "cg12934382",
                "cg08468401",
                "cg20816447",
                "cg00573770",
                "cg06335143",
                "cg06155229",
                "cg03032497",
                "cg06619077",
                "cg17804348",
                "cg00329615",
                "cg23479922",
                "cg10501210",
                "cg19991948",
                "cg27312979",
                "cg23186333",
                "cg25413977",
                "cg22078805",
                "cg17621438",
                "cg21878650",
                "cg04503319",
                "cg09809672",
            ],
            "default_imputation": "none",
        },
    },
    "GPAge71": {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00211-4",
        "output": "Age (Years)",
        "model": {
            "type": "GPAgeModel",
            "file": "GP-age_model_71_cpgs.json.zip",
            "sites": [
                "cg16867657",
                "cg22454769",
                "cg06639320",
                "cg04875128",
                "cg19283806",
                "cg24724428",
                "cg07553761",
                "cg24079702",
                "cg14556683",
                "cg07547549",
                "cg08128734",
                "cg23500537",
                "cg12934382",
                "cg08468401",
                "cg05404236",
                "cg20816447",
                "cg17110586",
                "cg24466241",
                "cg18473521",
                "cg00573770",
                "cg06335143",
                "cg06155229",
                "cg03032497",
                "cg12899747",
                "cg06619077",
                "cg17804348",
                "cg00329615",
                "cg23479922",
                "cg09017434",
                "cg10501210",
                "cg19991948",
                "cg03738025",
                "cg27312979",
                "cg14766700",
                "cg23186333",
                "cg21184711",
                "cg22730004",
                "cg19421125",
                "cg15894389",
                "cg10835286",
                "cg25413977",
                "cg11807280",
                "cg22078805",
                "cg02872426",
                "cg00303541",
                "cg04295144",
                "cg19729744",
                "cg24794228",
                "cg17621438",
                "cg05017994",
                "cg21878650",
                "cg15804973",
                "cg07797372",
                "cg27152890",
                "cg17471939",
                "cg18887458",
                "cg15243034",
                "cg04503319",
                "cg03350900",
                "cg09809672",
                "cg14577707",
                "cg24892069",
                "cg00602811",
                "cg05991454",
                "cg23126342",
                "cg09988805",
                "cg19344626",
                "cg20059012",
                "cg21826784",
                "cg07164639",
                "cg11693709",
            ],
            "default_imputation": "none",
        },
    },
    "GPAgeA": {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00211-4",
        "output": "Age (Years)",
        "model": {
            "type": "GPAgeModel",
            "file": "GP-age_model_a.json.zip",
            "sites": [
                "cg16867657",
                "cg19283806",
                "cg24724428",
                "cg07547549",
                "cg05404236",
                "cg06155229",
                "cg03738025",
                "cg14766700",
                "cg23186333",
                "cg21184711",
                "cg15894389",
                "cg02872426",
                "cg15804973",
                "cg18887458",
                "cg15243034",
                "cg03350900",
                "cg23126342",
                "cg07164639",
            ],
            "default_imputation": "none",
        },
    },
    "GPAgeB": {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00211-4",
        "output": "Age (Years)",
        "model": {
            "type": "GPAgeModel",
            "file": "GP-age_model_b.json.zip",
            "sites": [
                "cg22454769",
                "cg06639320",
                "cg07553761",
                "cg24079702",
                "cg08128734",
                "cg12934382",
                "cg08468401",
                "cg20816447",
                "cg24466241",
                "cg00573770",
                "cg06335143",
                "cg03032497",
                "cg12899747",
                "cg06619077",
                "cg17804348",
                "cg00329615",
                "cg10501210",
                "cg22730004",
                "cg10835286",
                "cg25413977",
                "cg11807280",
                "cg22078805",
                "cg00303541",
                "cg19729744",
                "cg07797372",
                "cg17471939",
                "cg09809672",
                "cg14577707",
                "cg00602811",
                "cg05991454",
                "cg09988805",
                "cg21826784",
            ],
            "default_imputation": "none",
        },
    },
    "GPAgeC": {
        "year": 2023,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00211-4",
        "output": "Age (Years)",
        "model": {
            "type": "GPAgeModel",
            "file": "GP-age_model_c.json.zip",
            "sites": [
                "cg04875128",
                "cg14556683",
                "cg23500537",
                "cg17110586",
                "cg18473521",
                "cg23479922",
                "cg09017434",
                "cg19991948",
                "cg27312979",
                "cg19421125",
                "cg04295144",
                "cg24794228",
                "cg17621438",
                "cg05017994",
                "cg21878650",
                "cg27152890",
                "cg04503319",
                "cg24892069",
                "cg19344626",
                "cg20059012",
                "cg11693709",
            ],
            "default_imputation": "none",
        },
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
    "DNAmClockCortical": {
        "year": 2020,
        "species": "Human",
        "tissue": "Human Cortex",
        "source": "https://doi.org/10.1093/brain/awaa334",
        "output": "Human Cortex Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "DNAmClockCortical.csv",
            "transform": lambda sum: anti_trafo(sum + 0.577682570446177),
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
    "VidalBralo": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.3389/fgene.2016.00126",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "VidalBralo.csv",
            "transform": lambda sum: sum + 84.7,
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
    "EpiTOC1": {
        "year": 2016,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/s13059-016-1064-3",
        "output": "Stem Cell Division Rate",
        "model": {
            "type": "LinearMethylationModel",
            "file": "EpiTOC1.csv",
        },
    },
    "EpiTOC2": {
        "year": 2020,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/s13073-020-00752-3",
        "output": "Stem Cell Division Rate",
        "model": {
            "type": "EpiTOC2Model",
            "file": "EpiTOC2.csv",
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
    "StocZ": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1038/s43587-024-00600-8",
        "output": "Mortality Risk",
        "model": {
            "type": "LinearMethylationModel",
            "file": "StocZ.csv",
            "transform": lambda sum: sum + 64.8077188694894,
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
    "StocP": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1038/s43587-024-00600-8",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "StocP.csv",
            "transform": lambda sum: sum + 92.8310813279039,
        },
    },
    "Mayne": {
        "year": 2016,
        "species": "Human",
        "tissue": "Placenta",
        "source": "https://doi.org/10.2217/epi-2016-0103",
        "output": "Gestational age",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Mayne.csv",
            "transform": lambda sum: sum + 24.99026,
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
    "AltumAge": {
        "year": 2022,
        "species": "Human",
        "tissue": "Multi-tissue",
        "source": "https://doi.org/10.1038/s41514-022-00085-y",
        "output": "Age (Years)",
        "model": {
            "type": "AltumAgeModel",
            "file": "AltumAge.csv",
            "weights": "AltumAge.pt",
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
    "StocH": {
        "year": 2024,
        "species": "Human",
        "tissue": "Multi-tissue",
        "source": "https://doi.org/10.1038/s43587-024-00600-8",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "StocH.csv",
            "transform": lambda sum: sum + 59.8015666314217,
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
    "Weidner": {
        "year": 2014,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1186/gb-2014-15-2-r24",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Weidner.csv",
            "transform": lambda sum: sum + 38.0,
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
    "Garagnani": {
        "year": 2012,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://pubmed.ncbi.nlm.nih.gov/23061750/",
        "output": "Age (Years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Garagnani.csv",
            "transform": lambda sum: sum,
        },
    },
    "OrganAgeChronological": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.medrxiv.org/content/10.1101/2024.04.08.24305469v1",
        "output": "Mortality Risk by Organ",
        "model": {
            "type": "LinearMultipartProteomicModel",
            "preprocess": olink_standardization_preprocess(
                "reference/olink3000_deviations.csv"
            ),
            "file": "OrganAgeChronological.csv",
            "default_imputation": "none",
        },
    },
    "OrganAgeMortality": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.medrxiv.org/content/10.1101/2024.04.08.24305469v1",
        "output": "Age (Years) by Organ",
        "model": {
            "type": "LinearMultipartProteomicModel",
            "preprocess": olink_standardization_preprocess(
                "reference/olink3000_deviations.csv"
            ),
            "file": "OrganAgeMortality.csv",
            "default_imputation": "none",
        },
    },
    "OrganAge1500Chronological": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.medrxiv.org/content/10.1101/2024.04.08.24305469v1",
        "output": "Mortality Risk by Organ",
        "model": {
            "type": "LinearMultipartProteomicModel",
            "preprocess": olink_standardization_preprocess(
                "reference/olink1500_deviations.csv"
            ),
            "file": "OrganAge1500Chronological.csv",
            "default_imputation": "none",
        },
    },
    "OrganAge1500Mortality": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://www.medrxiv.org/content/10.1101/2024.04.08.24305469v1",
        "output": "Age (Years) by Organ",
        "model": {
            "type": "LinearMultipartProteomicModel",
            "preprocess": olink_standardization_preprocess(
                "reference/olink1500_deviations.csv"
            ),
            "file": "OrganAge1500Mortality.csv",
            "default_imputation": "none",
        },
    },
    "Bohlin": {
        "year": 2017,
        "species": "Human",
        "tissue": "Cord Blood",
        "source": "https://doi.org/10.1186/s13059-016-1063-4",
        "output": "Age (days)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Bohlin.csv",
            "transform": lambda sum: sum + 277.2421,
        },
    },
    "Bocklandt": {
        "year": 2011,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1371/journal.pone.0014821",
        "output": "Age (years)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Bocklandt.csv",
        },
    },
    "MiAge": {
        "year": 2018,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://doi.org/10.1080/15592294.2017.1389361",
        "output": "Mitotic Age (Cell Divisions)",
        "model": {
            "type": "MiAgeModel",
            "file": "MiAge.csv",
        },
    },
    "Bohlin": {
        "year": 2017,
        "species": "Human",
        "tissue": "Cord Blood",
        "source": "https://doi.org/10.1186/s13059-016-1063-4",
        "output": "Age (days)",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Bohlin.csv",
            "transform": lambda sum: sum + 277.2421,
        },
    },
    "PCHorvath1": {
        "year": 2022,
        "species": "Human",
        "tissue": "Multi-tissue",
        "source": "https://doi.org/10.1038/s43587-022-00248-2",
        "output": "Age (Years)",
        "model": {
            "type": "PCLinearTransformationModel",
            "file": "https://storage.googleapis.com/biolearn/PCClock/PCHorvath1_model.csv",
            "rotation": "https://storage.googleapis.com/biolearn/PCClock/PCHorvath1_rotation.csv",
            "center": "https://storage.googleapis.com/biolearn/PCClock/PCHorvath1_center.csv",
            "transform": lambda sum: anti_trafo(sum + 1.15834584357227),
            "default_imputation": "sesame_450k",
        },
    },
    "HurdleInflammAge": {
        "year": 2024,
        "species": "Human",
        "tissue": "Blood",
        "source": "https://hurdle.bio/our-science/",
        "output": "InflammAge",
        "model": {
            "type": "HurdleAPIModel",
            "use_production": True,
            "default_imputation": "none",
        },
        "usage": {
            "commercial": "Non-Commercial Use Only. Any commercial use, distribution, or resale of this information is strictly prohibited without prior written agreement from Chronomics Limited. See https://hurdle.bio/ for commercial information.",
            "non-commercial": "Register at https://dashboard.hurdle.bio/register then navigate to Developer > API keys",
        },
    },
    "Pasta": {
        "year": 2025,
        "species": "Human",
        "tissue": "Multi-tissue",
        "source": "https://www.biorxiv.org/content/10.1101/2025.06.04.657785v1.full",
        "output": "Age score",
        "model": {
            "type": "LinearTranscriptomicModel",
            "file": "Pasta.csv",
            "preprocess": preprocess_pasta,
            "transform": lambda sum: sum * -4.76348378687217
            - 0.0502893445253186 * -4.76348378687217,
        },
        "usage": {
            "commercial": "Free to use",
            "non-commercial": "Free to use",
        },
    },
    "REG": {
        "year": 2025,
        "species": "Human",
        "tissue": "Multi-tissue",
        "source": "https://www.biorxiv.org/content/10.1101/2025.06.04.657785v1.full",
        "output": "Age (years)",
        "model": {
            "type": "LinearTranscriptomicModel",
            "file": "REG.csv",
            "preprocess": preprocess_pasta,
            "transform": lambda sum: sum + 140.272578432562,
        },
        "usage": {
            "commercial": "Free to use",
            "non-commercial": "Free to use",
        },
    },
}


class AltumAgeNeuralNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        self.linear1 = nn.Linear(20318, 32)
        self.linear2 = nn.Linear(32, 32)
        self.linear3 = nn.Linear(32, 32)
        self.linear4 = nn.Linear(32, 32)
        self.linear5 = nn.Linear(32, 32)
        self.linear6 = nn.Linear(32, 1)

        self.bn1 = nn.BatchNorm1d(20318, eps=0.001, momentum=0.99)
        self.bn2 = nn.BatchNorm1d(32, eps=0.001, momentum=0.99)
        self.bn3 = nn.BatchNorm1d(32, eps=0.001, momentum=0.99)
        self.bn4 = nn.BatchNorm1d(32, eps=0.001, momentum=0.99)
        self.bn5 = nn.BatchNorm1d(32, eps=0.001, momentum=0.99)
        self.bn6 = nn.BatchNorm1d(32, eps=0.001, momentum=0.99)

    def forward(self, x):
        x = self.bn1(x)
        x = self.linear1(x)
        x = F.selu(x)

        x = self.bn2(x)
        x = self.linear2(x)
        x = F.selu(x)

        x = self.bn3(x)
        x = self.linear3(x)
        x = F.selu(x)

        x = self.bn4(x)
        x = self.linear4(x)
        x = F.selu(x)

        x = self.bn5(x)
        x = self.linear5(x)
        x = F.selu(x)

        x = self.bn6(x)
        x = self.linear6(x)
        return x


class AltumAgeModel:
    def __init__(self, weights_path=None, preprocess_file_path=None):
        self.model = AltumAgeNeuralNetwork()
        self.model.eval()
        self._weights_loaded = False
        self.center = None
        self.scale = None
        self.reference = None

        if preprocess_file_path:
            if not os.path.isabs(preprocess_file_path):
                preprocess_file_path = get_data_file(preprocess_file_path)
            preprocess_df = pd.read_csv(
                preprocess_file_path, index_col=0
            )  # Ensure CpG_site is the index
            self.center = torch.tensor(
                preprocess_df["center"].values, dtype=torch.float32
            )
            self.scale = torch.tensor(
                preprocess_df["scale"].values, dtype=torch.float32
            )
            self.reference = (
                preprocess_df.index.tolist()
            )  # Convert index to a list of CpG site names

        if weights_path:
            self.load_pytorch_weights(weights_path)

    def load_pytorch_weights(self, weights_path: str):
        p = weights_path
        if not os.path.isabs(p):
            p = get_data_file(p)
        if not os.path.exists(p):
            print(
                f"Warning: Weights file {p} does not exist. Model will be initialized with random weights."
            )
            self._weights_loaded = False
            return
        sd = torch.load(p, map_location="cpu", weights_only=False)
        if isinstance(sd, dict) and "state_dict" in sd:
            sd = sd["state_dict"]
        elif (
            isinstance(sd, dict)
            and "model" in sd
            and isinstance(sd["model"], dict)
        ):
            sd = sd["model"]
        self.model.load_state_dict(sd, strict=True)
        self._weights_loaded = True

    def predict(self, geo_data):
        """
        Predicts the output using the AltumAgeModel.

        Args:
            geo_data (GeoData): A GeoData object containing metadata and dnam attributes.

        Returns:
            pd.DataFrame: Predictions from the model with a column named "Predicted".
        """
        if self.reference is None or self.center is None or self.scale is None:
            raise RuntimeError("AltumAge preprocessing parameters missing.")
        if not self._weights_loaded:
            raise RuntimeError(
                "AltumAge used without loaded weights. Please provide AltumAge.pt file."
            )
        # Access the DNA methylation data from the GeoData object
        df = geo_data.dnam.fillna(0)
        # Ensure the input DataFrame contains all required CpG sites
        required_cpgs = self.methylation_sites()
        missing_cpgs = set(required_cpgs) - set(df.index)
        if missing_cpgs:
            print(
                f"{len(missing_cpgs)} CpGs missing in input; filling with 0."
            )

        # Align input data with the reference CpG sites
        df = df.reindex(
            self.reference, fill_value=0
        )  # Fill missing CpG sites with 0

        # # Convert input DataFrame to a PyTorch tensor
        X = torch.tensor(
            df.values.T, dtype=torch.float32
        )  # Transpose to match sample orientation

        if self.center is not None and self.scale is not None:
            # Scale the input using the center and scale values
            X = (X - self.center) / (self.scale + 1e-8)
        # Pass the preprocessed input through the neural network
        with torch.no_grad():
            preds = self.model(X).squeeze().numpy()

        # Return predictions as a DataFrame
        return pd.DataFrame(preds, index=df.columns, columns=["Predicted"])

    @classmethod
    def from_definition(cls, clock_definition):
        """
        Creates an instance of AltumAgeModel from a model definition.
        Args:
            clock_definition (dict): The model definition containing file paths.
        Returns:
            AltumAgeModel: An instance of the AltumAgeModel class.
        """
        model_def = clock_definition["model"]

        # Extract the file path for preprocessing
        preprocess_file = get_data_file(model_def["file"])
        weights_path = None

        try:
            weights_path = get_data_file(model_def["weights"])
        except Exception:
            weights_path = None
        return cls(
            weights_path=weights_path, preprocess_file_path=preprocess_file
        )

    def methylation_sites(self):
        return list(self.reference)


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
        cell_prop_df = pd.DataFrame(
            cell_prop, columns=sample_names, index=cell_types
        )

        # Return samples as rows to match other model outputs
        return cell_prop_df.T

    # returns required methylation sites
    def methylation_sites(self):
        return list(self.reference.index)


class LinearModel:
    def __init__(
        self,
        coefficient_file_or_df,
        transform,
        preprocess=None,
        **details,
    ) -> None:
        self.transform = transform
        self.preprocess = preprocess
        self.details = details

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
            center_file=model_def.get("center"),
            rotation_file=model_def.get("rotation"),
            **{k: v for k, v in clock_definition.items() if k != "model"},
        )

    def predict(self, geo_data):
        matrix_data = self._get_data_matrix(geo_data)
        matrix_data = self.preprocess(matrix_data)
        self._validate_required_features(matrix_data)
        matrix_data.loc["intercept"] = 1

        # Join the coefficients and dnam_data on the index
        model_df = self.coefficients.join(matrix_data, how="inner")

        # Vectorized multiplication: multiply CoefficientTraining with all columns of dnam_data
        coefficient_column = (
            "CoefficientTraining"
            if "CoefficientTraining" in model_df.columns
            else "Weight"
        )

        result = (
            model_df.iloc[:, 1:]
            .multiply(model_df[coefficient_column], axis=0)
            .sum(axis=0)
        )

        # Return as a DataFrame
        return result.apply(self.transform).to_frame(name="Predicted")

    def _validate_required_features(self, matrix_data):
        return

    def _get_data_matrix(self, geo_data):
        raise NotImplementedError()


class LinearMethylationModel(LinearModel):
    _MISSING_CPG_PREVIEW_LIMIT = 5

    def _get_data_matrix(self, geo_data):
        return geo_data.dnam

    def _validate_required_features(self, matrix_data):
        required_cpgs = self.methylation_sites()
        missing_cpgs = sorted(set(required_cpgs) - set(matrix_data.index))
        if not missing_cpgs:
            return

        model_name = self.details.get("name")
        model_label = f" for model '{model_name}'" if model_name else ""
        preview_limit = self._MISSING_CPG_PREVIEW_LIMIT
        preview = ", ".join(missing_cpgs[:preview_limit])
        remaining = len(missing_cpgs) - preview_limit
        if remaining > 0:
            preview = (
                f"showing first {preview_limit}: {preview} (+{remaining} more)"
            )

        raise ValueError(
            "Missing required CpG sites"
            f"{model_label} ({len(missing_cpgs)}/{len(required_cpgs)}): "
            f"{preview}. "
            "Provide methylation data with these CpGs or use an imputation method that includes them."
        )

    def methylation_sites(self):
        unique_vars = set(self.coefficients.index) - {"intercept"}
        return list(unique_vars)


class PCLinearTransformationModel(LinearModel):
    """
    Transforms methylation data into principal component (PC) space using
    provided center and rotation files before performing linear regression.
    """

    def __init__(
        self,
        coefficient_file_or_df,
        transform,
        preprocess=None,
        center_file=None,
        rotation_file=None,
        **details,
    ):
        """
        Initialize the PCLinearTransformationModel.

        Args:
            coefficient_file_or_df: Path to coefficient file or DataFrame
            transform: Transformation function to apply to predictions
            preprocess: Optional preprocessing function
            center_file (str): Path to the CSV file containing CpG means for centering.
            rotation_file (str): Path to the CSV file containing PCA loadings.
            **details: Additional details passed to the parent class.
        """
        # Pass center_file and rotation_file as part of **details
        details["center"] = center_file
        details["rotation"] = rotation_file
        super().__init__(
            coefficient_file_or_df, transform, preprocess, **details
        )
        self.center_file = center_file or details.get("center")
        self.rotation_file = rotation_file or details.get("rotation")
        self.center_ = None
        self.rotation = None

    def _load_data(self):
        """
        Load the center and rotation data from the provided files.
        """
        self.center_ = pd.read_csv(
            get_data_file(self.center_file), index_col=0
        ).iloc[:, 0]
        # Load rotation data if not already loaded
        self.rotation = pd.read_csv(
            get_data_file(self.rotation_file), index_col=None
        )
        self.rotation.set_index(self.center_.index, inplace=True)

        # Align indices to ensure consistency
        common_indices = self.center_.index.intersection(self.rotation.index)
        if len(common_indices) == 0:
            raise ValueError(
                "No common CpG sites found between center and rotation data."
            )
        self.center_ = self.center_.loc[common_indices]
        self.rotation = self.rotation.loc[common_indices]

    def _get_data_matrix(self, geo_data):
        """
        Apply PCA transformation to the DNA methylation data.

        Args:
            geo_data (GeoData): The input GeoData object containing methylation data.

        Returns:
            pd.DataFrame: The PC-transformed data matrix.
        """
        # Load the center and rotation data
        self._load_data()

        meth = geo_data.dnam.copy()

        # Intersect with CpGs in reference files
        intersecting_cpgs = meth.index.intersection(self.center_.index)

        meth = meth.loc[intersecting_cpgs]
        center_ = self.center_.loc[intersecting_cpgs]
        rotation = self.rotation.loc[intersecting_cpgs]

        # Fill NaN values with center values (so they become 0 after centering)
        meth = meth.T.fillna(center_).T

        # Center the data and apply PCA transformation
        X_centered = meth.subtract(center_, axis=0)
        PCs = X_centered.T.dot(rotation)  # (samples × PCs)
        return PCs.T  # (PCs × samples)

    def methylation_sites(self):
        """
        Return the list of required CpG sites.
        """
        # Load the center data if not already loaded
        if self.center_ is None:
            if not self.center_file:
                raise ValueError("Center file path is not provided.")
            self.center_ = pd.read_csv(
                get_data_file(self.center_file), index_col=0
            ).iloc[:, 0]
        return list(self.center_.index)


class LinearTranscriptomicModel(LinearModel):
    def _get_data_matrix(self, geo_data):
        return geo_data.rna


class GrimageModel:
    def __init__(self, coefficient_file, **details):
        self.coefficients = pd.read_csv(
            get_data_file(coefficient_file), index_col=0
        )
        self.details = details

    @classmethod
    def from_definition(cls, clock_definition):
        model_def = clock_definition["model"]
        return cls(
            model_def["file"],
            **{k: v for k, v in clock_definition.items() if k != "model"},
        )

    def predict(self, geo_data):
        if "sex" not in geo_data.metadata:
            raise ValueError(
                "GrimAge requires 'sex' column in metadata. "
                "If sex is unknown, you can predict it using the SexEstimation model:\n"
                "  from biolearn.model_gallery import ModelGallery\n"
                "  gallery = ModelGallery()\n"
                "  sex_pred = gallery.get('SexEstimation').predict(your_data)\n"
                "  # Then add sex to metadata based on predictions"
            )

        if "age" not in geo_data.metadata:
            raise ValueError(
                "GrimAge requires 'age' column in metadata (numeric age in years)"
            )

        # Check for NaN sex values
        sex_values = geo_data.metadata["sex"]
        if sex_values.isna().any():
            nan_samples = sex_values[sex_values.isna()].index.tolist()
            raise ValueError(
                f"GrimAge cannot process samples with unknown sex: {nan_samples[:3]}.\n"
                f"Either exclude these samples or predict sex using SexEstimation model."
            )

        # Check for sample mismatch between methylation matrix and metadata
        dnam_samples = set(geo_data.dnam.columns)
        metadata_samples = set(geo_data.metadata.index)
        missing_metadata = dnam_samples - metadata_samples

        if missing_metadata:
            raise ValueError(
                f"Methylation data contains samples without metadata: {list(missing_metadata)[:3]}.\n"
                f"Ensure all samples in methylation matrix have corresponding metadata entries."
            )

        df = geo_data.dnam

        # Transposing metadata so that its structure aligns with dnam (columns as samples)
        transposed_metadata = geo_data.metadata.transpose()

        # Add metadata rows to dnam DataFrame
        df.loc["Age"] = transposed_metadata.loc["age"]
        df.loc["Female"] = transposed_metadata.loc["sex"].apply(
            lambda x: 1 if x == 0 else 0
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
            lambda x: 1 if x == 0 else 0
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

        # Check for NaN results which indicate missing data
        if result.isna().any():
            nan_samples = result[result.isna()].index.tolist()
            raise ValueError(
                f"Missing methylation data for required CpG sites in samples: {nan_samples[:3]}. "
                f"Ensure all required methylation sites are present in your data."
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


class LinearMultipartProteomicModel:
    def __init__(
        self, coefficient_file_or_df, transform, preprocess, **details
    ):
        self.details = details
        self.transform = transform
        self.preprocess = preprocess
        if isinstance(coefficient_file_or_df, pd.DataFrame):
            self.coefficients = coefficient_file_or_df
        else:
            file_path = (
                coefficient_file_or_df
                if os.path.isabs(coefficient_file_or_df)
                else get_data_file(coefficient_file_or_df)
            )
            self.coefficients = pd.read_csv(file_path)

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
        mat = geo_data.protein_olink
        if mat is None or mat.empty:
            raise ValueError(
                "No olink proteomic data provided: 'geo_data.protein_olink' is None or empty."
            )

        # Apply preprocessing
        mat = self.preprocess(mat)
        results = {}
        for tissue, grp in self.coefficients.groupby("Tissue"):
            # intercept
            intercepts = grp.loc[
                grp["Protein"].str.lower() == "intercept", "Coefficient"
            ]
            intercept = (
                float(intercepts.iloc[0]) if not intercepts.empty else 0.0
            )

            # protein coefficients
            coef_ser = grp.loc[
                grp["Protein"].str.lower() != "intercept"
            ].set_index("Protein")["Coefficient"]

            # align proteins as columns, fill missing
            aligned = mat.reindex(columns=coef_ser.index).fillna(0)

            # compute prediction per sample
            pred = aligned.dot(coef_ser) + intercept
            results[tissue] = pred

        # Apply transformation to results
        return self.transform(pd.DataFrame(results))

    def methylation_sites(self):
        return []


class SexEstimationModel:
    def __init__(self, coeffecient_file, **details):
        self.coefficients = pd.read_csv(
            get_data_file(coeffecient_file), index_col=0, low_memory=False
        )
        self.details = details

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


class EpiTOC2Model:
    def __init__(self, reference_file):
        self.reference = pd.read_csv(
            get_data_file(reference_file), index_col=0
        )
        self.delta = self.reference["delta"].to_numpy().reshape(-1, 1)
        self.beta0 = self.reference["beta0"].to_numpy().reshape(-1, 1)
        self.CpG_names = self.reference.index.to_numpy()

    @classmethod
    def from_definition(cls, clock_definition):
        model_def = clock_definition["model"]
        return cls(model_def["file"])

    def predict(self, geo_data):
        dnam = geo_data.dnam

        map_idx = dnam.index.get_indexer(self.CpG_names)
        rep_mask = map_idx != -1

        rows = map_idx[rep_mask]
        beta = dnam.values[rows, :].astype(float)
        delta = self.delta[rep_mask, :]
        beta0 = self.beta0[rep_mask, :]
        k = float(rows.size)

        beta = np.where(np.isnan(beta), 0.0, beta)

        with np.errstate(divide="ignore", invalid="ignore"):
            contrib = (beta - beta0) / (delta * (1.0 - beta0))

        vals = 2.0 * (np.sum(contrib, axis=0) / k)

        return pd.DataFrame(vals, index=dnam.columns, columns=["Predicted"])

    def methylation_sites(self):
        return list(self.CpG_names)


class HurdleAPIModel:
    """
    Hurdle Bio InflammAge model via API.

    Non-Commercial Use Only:
    This inflammaging calculation is provided for personal and educational purposes only.
    Any commercial use, distribution, or resale of this information is strictly prohibited
    without prior written agreement from Chronomics Limited. See https://hurdle.bio/ for
    further commercial information.

    Requires API credentials from https://dashboard.hurdle.bio/register
    Navigate to Developer > API keys after registration.
    """

    # API Configuration
    SANDBOX_BASE_URL = "https://api.sandbox.hurdle.bio"
    PRODUCTION_BASE_URL = "https://api.hurdle.bio"
    API_ENDPOINT_PATH = "/predict/v1/"
    DEFAULT_TIMEOUT = 30  # seconds
    MAX_RETRIES = 3
    RETRY_DELAY = 1  # seconds

    def __init__(
        self,
        api_key: Optional[str] = None,
        use_production: bool = False,
        timeout: Optional[int] = None,
        **details: Any,
    ) -> None:
        self.details = details
        self.api_key = api_key or os.environ.get("HURDLE_API_KEY")
        self.timeout = timeout or self.DEFAULT_TIMEOUT
        self._consent_given = False

        if not self.api_key:
            raise ValueError(
                "API key required. Either pass api_key parameter or set HURDLE_API_KEY environment variable.\n"
                "Get your API key at: https://dashboard.hurdle.bio/register"
            )

        base_url = (
            self.PRODUCTION_BASE_URL
            if use_production
            else self.SANDBOX_BASE_URL
        )
        self.api_endpoint = f"{base_url}{self.API_ENDPOINT_PATH}"

        # Load required CpG sites
        try:
            cpg_file = get_data_file("Hurdle_CpGs.csv")
            self.required_cpgs = pd.read_csv(cpg_file, encoding="ISO-8859-1")[
                "ProbeID"
            ].tolist()
        except (FileNotFoundError, KeyError, pd.errors.EmptyDataError) as e:
            self.required_cpgs = None
            warnings.warn(
                f"Could not load Hurdle CpG sites: {e}. Will use all available CpG sites."
            )

    @classmethod
    def from_definition(
        cls, clock_definition: Dict[str, Any]
    ) -> "HurdleAPIModel":
        model_def = clock_definition["model"]
        return cls(
            use_production=model_def.get("use_production", False),
            **{k: v for k, v in clock_definition.items() if k != "model"},
        )

    def _get_consent(self, require_consent: bool) -> None:
        if require_consent and not self._consent_given:
            print(
                "\nWARNING: This will send methylation data to Hurdle Bio's servers."
            )
            print(
                "Please ensure you have permission to share this data with third parties."
            )
            response = (
                input("Do you consent to send this data? (yes/no): ")
                .lower()
                .strip()
            )

            if response != "yes":
                raise ValueError(
                    "User consent required to send data to external API"
                )
            self._consent_given = True

    def predict(
        self,
        data: Union[pd.DataFrame, "GeoData"],
        require_consent: bool = True,
    ) -> pd.Series:
        self._get_consent(require_consent)

        # Extract methylation data
        if isinstance(data, GeoData):
            dnam = data.dnam
            metadata = data.metadata
        else:
            dnam = data
            metadata = pd.DataFrame(index=dnam.columns)

        # Filter and impute CpG sites
        if self.required_cpgs:
            missing_cpgs = set(self.required_cpgs) - set(dnam.index)
            available_cpgs = [
                cpg for cpg in self.required_cpgs if cpg in dnam.index
            ]

            if not available_cpgs:
                raise ValueError("No required CpG sites found in data")

            # Check for missing CpGs and raise error if any are missing
            if missing_cpgs:
                n_missing = len(missing_cpgs)
                n_required = len(self.required_cpgs)
                sample_missing = list(missing_cpgs)[:5]

                raise ValueError(
                    f"Missing {n_missing}/{n_required} required CpG sites. "
                    f"Examples: {sample_missing}. "
                    f"Please impute missing values before calling predict() or use "
                    f"ModelGallery.get('HurdleInflammAge', imputation_method='averaging')"
                )

            filtered_dnam = dnam.loc[available_cpgs]
        else:
            filtered_dnam = dnam

        # Round values
        filtered_dnam = filtered_dnam.round(3)

        try:
            # Get presigned URL
            headers = {"x-api-key": self.api_key}
            response = requests.post(
                f"{self.api_endpoint}upload_link",
                json={"fileType": "beta_matrix"},
                headers=headers,
                timeout=self.timeout,
            )

            if response.status_code != 200:
                raise Exception(
                    f"Failed to get upload URL: {response.status_code} - {response.text}"
                )

            upload_info = response.json()

            # Validate response structure
            if not isinstance(upload_info, dict):
                raise ValueError(
                    f"Invalid API response format: expected dict, got {type(upload_info)}"
                )
            if (
                "uploadUrl" not in upload_info
                or "requestId" not in upload_info
            ):
                raise ValueError(
                    f"Missing required fields in API response: {upload_info.keys()}"
                )

            # Upload beta matrix
            csv_buffer = io.StringIO()
            filtered_dnam.to_csv(csv_buffer)
            csv_buffer.seek(0)

            upload_response = requests.put(
                upload_info["uploadUrl"],
                data=csv_buffer.getvalue().encode("utf-8"),
                timeout=self.timeout,
            )

            if upload_response.status_code != 200:
                raise Exception(
                    f"Failed to upload data: {upload_response.status_code}"
                )

            # Prepare metadata
            meta_list = []
            for sample_id in filtered_dnam.columns:
                sample_meta = {"barcode": sample_id, "tissue": "whole_blood"}

                if sample_id in metadata.index:
                    row = metadata.loc[sample_id]
                    if "age" in row:
                        sample_meta["chronologicalAge"] = float(row["age"])
                    if "sex" in row:
                        sex_value = row["sex"]
                        # Biolearn standard: 0=female, 1=male (metadata_standard.rst)
                        # Strict validation - only accept standard values
                        if sex_value == 0:
                            sample_meta["sex"] = "f"
                        elif sex_value == 1:
                            sample_meta["sex"] = "m"
                        else:
                            raise ValueError(
                                f"Invalid sex value '{sex_value}' for sample "
                                f"'{sample_id}'. Biolearn standard requires: "
                                f"0 (female) or 1 (male). "
                                f"See metadata_standard.rst for details."
                            )

                meta_list.append(sample_meta)

            # Get predictions
            prediction_request = {
                "requestId": upload_info["requestId"],
                "biomarker": "inflammage",
                "arrayType": "hurdle",
                "metadata": meta_list,
            }

            pred_response = requests.post(
                self.api_endpoint, json=prediction_request, headers=headers
            )

            if pred_response.status_code != 200:
                raise Exception(
                    f"Prediction failed: {pred_response.status_code} - {pred_response.text}"
                )

            # Extract predictions
            results = pred_response.json()

            # Validate response structure
            if not isinstance(results, dict):
                raise ValueError(
                    f"Invalid API response format: expected dict, got {type(results)}"
                )
            if "data" not in results:
                raise ValueError(
                    f"Missing 'data' field in API response: {results.keys()}"
                )
            if not isinstance(results["data"], list):
                raise ValueError(
                    f"Invalid 'data' field: expected list, got {type(results['data'])}"
                )

            predictions = {}
            for item in results.get("data", []):
                if (
                    not isinstance(item, dict)
                    or "barcode" not in item
                    or "prediction" not in item
                ):
                    raise ValueError(f"Invalid prediction item format: {item}")
                predictions[item["barcode"]] = float(item["prediction"])

            return pd.Series(predictions, name="inflammage")

        except requests.exceptions.RequestException as e:
            raise Exception(f"Network error: {str(e)}")
        except Exception as e:
            raise Exception(f"API error: {str(e)}")

    def methylation_sites(self):
        """Return list of required CpG sites for imputation compatibility."""
        return self.required_cpgs if self.required_cpgs else []


class GPAgeModel:
    """Gaussian Process regression model for age prediction (GP-age clock)."""

    def __init__(self, model_file, sites):
        self.model_file = model_file
        self._sites = sites
        self._model = None
        self._training_means = None

    def _load_model(self):
        if self._model is not None:
            return
        try:
            import GPy
        except ImportError:
            raise ImportError(
                "GPy is required for GP-age models. "
                "Install with: pip install biolearn[gpage]"
            )
        model_path = get_data_file(self.model_file)
        self._model = GPy.models.GPRegression.load_model(model_path)
        self._training_means = np.mean(self._model.X, axis=0)

    @classmethod
    def from_definition(cls, clock_definition):
        model_def = clock_definition["model"]
        return cls(model_def["file"], model_def["sites"])

    def predict(self, geo_data):
        self._load_model()
        dnam = geo_data.dnam
        methylation_data = pd.DataFrame(
            np.nan, index=dnam.columns, columns=self._sites
        )
        existing = [s for s in self._sites if s in dnam.index]
        for site in existing:
            methylation_data[site] = dnam.loc[site].values
        X = methylation_data.values
        if np.isnan(X).any():
            for i, mean_val in enumerate(self._training_means):
                X[np.isnan(X[:, i]), i] = mean_val
        predictions = self._model.predict(X)[0].squeeze()
        return pd.DataFrame(
            predictions, index=dnam.columns, columns=["Predicted"]
        )

    def methylation_sites(self):
        return self._sites


class ImputationDecorator:
    def __init__(self, clock, imputation_method):
        self.clock = clock
        self.imputation_method = imputation_method

    def predict(self, geo_data):
        needed_cpgs = self.clock.methylation_sites()
        dnam_data_imputed = self.imputation_method(geo_data.dnam, needed_cpgs)

        geo_copy = geo_data.copy()
        geo_copy.dnam = dnam_data_imputed
        return self.clock.predict(geo_copy)

    # Forwarding other methods and attributes to the clock
    def __getattr__(self, name):
        return getattr(self.clock, name)


class MiAgeModel:
    """
    MiAge (Mitotic Age) clock implementation.

    Based on Youn & Wang (2018): "The MiAge Calculator: a DNA methylation-based
    mitotic age calculator of human tissue types."

    MiAge estimates the number of cell divisions (mitotic age) using an optimization
    approach with site-specific parameters (b, c, d) for 268 CpG sites.

    Formula: For each sample, find n (mitotic age) that minimizes:
        sum((c + b^(n-1) * d - beta)^2)
    where beta is the observed methylation value for each CpG site.
    """

    def __init__(self, parameter_file):
        """
        Initialize MiAge model with site-specific parameters.

        Parameters
        ----------
        parameter_file : str
            Path to CSV file containing CpG, b, c, d columns
        """
        params = pd.read_csv(get_data_file(parameter_file))
        self.params = params.set_index("CpG")
        self.b_params = self.params["b"]
        self.c_params = self.params["c"]
        self.d_params = self.params["d"]
        self.cpg_sites = self.params.index.tolist()

    @classmethod
    def from_definition(cls, clock_definition):
        """Create MiAge model from clock definition"""
        model_def = clock_definition["model"]
        return cls(model_def["file"])

    def _optimize_sample(self, beta_values):
        """
        Optimize mitotic age for a single sample using scipy.optimize.

        Parameters
        ----------
        beta_values : numpy.ndarray
            Methylation beta values for CpG sites

        Returns
        -------
        float
            Estimated mitotic age (number of cell divisions)
        """
        # Filter to valid (non-NaN) values
        valid_mask = ~np.isnan(beta_values)
        beta_valid = beta_values[valid_mask]
        b_valid = self.b_params.values[valid_mask]
        c_valid = self.c_params.values[valid_mask]
        d_valid = self.d_params.values[valid_mask]

        def objective(n):
            """MiAge objective function: sum((c + b^(n-1) * d - beta)^2)"""
            return np.sum(
                (c_valid + b_valid ** (n - 1) * d_valid - beta_valid) ** 2
            )

        # Use multiple starting points for robustness
        best_result = None
        best_fun = np.inf

        for start in [100, 500, 1000, 2000]:
            result = minimize_scalar(
                objective,
                bounds=(10, 10000),
                method="bounded",
                options={"xatol": 1e-8},
            )
            if result.fun < best_fun:
                best_fun = result.fun
                best_result = result

        return best_result.x

    def predict(self, geo_data):
        """
        Predict mitotic age for methylation samples.

        Parameters
        ----------
        geo_data : GeoData
            Object containing methylation data

        Returns
        -------
        pd.DataFrame
            Predicted mitotic ages for each sample
        """
        methylation_data = geo_data.dnam

        # Get shared CpG sites
        shared_sites = methylation_data.index.intersection(self.cpg_sites)

        if len(shared_sites) == 0:
            raise ValueError("No overlapping CpG sites found for MiAge clock")

        if len(shared_sites) < len(self.cpg_sites) * 0.5:
            print(
                f"WARNING: Only {len(shared_sites)}/{len(self.cpg_sites)} "
                f"CpG sites available ({100*len(shared_sites)/len(self.cpg_sites):.1f}%)"
            )

        # Align data with model CpGs (keeping only shared sites)
        # Need to reindex params to match shared sites order
        shared_params = self.params.loc[shared_sites]
        self.b_params = shared_params["b"]
        self.c_params = shared_params["c"]
        self.d_params = shared_params["d"]

        X = methylation_data.loc[shared_sites].values

        # Calculate MiAge for each sample
        predictions = []
        for i in range(X.shape[1]):
            beta_sample = X[:, i]
            miage = self._optimize_sample(beta_sample)
            predictions.append(miage)

        return pd.DataFrame(
            predictions, index=methylation_data.columns, columns=["Predicted"]
        )

    def methylation_sites(self):
        """Return list of required CpG sites"""
        return self.cpg_sites


def single_sample_clock(clock_function, data):
    return clock_function(data).iloc[0, 0]
