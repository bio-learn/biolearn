import pandas as pd
import numpy as np


def phenotypic_age(df):
    """Calculate phenotypic age from blood biomarkers.

    Based on: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940111/
    """
    pheno_coefs = {
        "age": 0.0804,
        "albumin": -0.034,
        "creatinine": 0.0095,
        "glucose": 0.1953,
        "c_reactive_protein": 0.0954,
        "lymphocyte_percent": -0.012,
        "mean_cell_volume": 0.0268,
        "red_blood_cell_distribution_width": 0.3356,
        "alkaline_phosphate": 0.00188,
        "white_blood_cell_count": 0.0554,
    }

    constant = -19.9067
    gamma = 0.0077
    cs = [141.50225, -0.00553, 0.090165]

    # Vectorized calculation - no DataFrame modifications
    pheno = (
        df["age"] * 0.0804
        + df["albumin"] * -0.034
        + df["creatinine"] * 0.0095
        + df["glucose"] * 0.1953
        + np.log(df["c_reactive_protein"]) * 0.0954
        + df["lymphocyte_percent"] * -0.012
        + df["mean_cell_volume"] * 0.0268
        + df["red_blood_cell_distribution_width"] * 0.3356
        + df["alkaline_phosphate"] * 0.00188
        + df["white_blood_cell_count"] * 0.0554
    )

    mortality_score = 1 - np.exp(
        (-np.exp(constant + pheno) * ((np.exp(120 * gamma) - 1) / gamma))
    )

    phenotypic_age = (
        cs[0] + np.log(cs[1] * np.log(1.00001 - mortality_score)) / cs[2]
    )

    return phenotypic_age
