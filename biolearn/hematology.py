import pandas as pd
import numpy as np


def phenotypic_age(df):
    """Calculate phenotypic age from blood biomarkers.

    Based on: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940111/

    Coefficients follow Levine et al. (2018), Aging 10(4):573-591,
    Table 1 / Supplementary Table S1, with gamma = 0.0076927 from
    Supplement 1. The alkaline phosphatase weight (0.00188) is the unrounded
    value printed in the correction to Liu et al. (2018), PLoS Med
    16(2):e1002760 (Table 1 rounds it to 0.0019).
    """
    pheno_coefs = {
        "age": 0.0804,
        "albumin": -0.0336,
        "creatinine": 0.0095,
        "glucose": 0.1953,
        "c_reactive_protein": 0.0954,
        "lymphocyte_percent": -0.012,
        "mean_cell_volume": 0.0268,
        "red_blood_cell_distribution_width": 0.3306,
        "alkaline_phosphate": 0.00188,
        "white_blood_cell_count": 0.0554,
    }

    constant = -19.9067
    gamma = 0.0076927
    cs = [141.50225, -0.00553, 0.090165]

    # Vectorized calculation - no DataFrame modifications
    pheno = (
        df["age"] * 0.0804
        + df["albumin"] * -0.0336
        + df["creatinine"] * 0.0095
        + df["glucose"] * 0.1953
        + np.log(df["c_reactive_protein"]) * 0.0954
        + df["lymphocyte_percent"] * -0.012
        + df["mean_cell_volume"] * 0.0268
        + df["red_blood_cell_distribution_width"] * 0.3306
        + df["alkaline_phosphate"] * 0.00188
        + df["white_blood_cell_count"] * 0.0554
    )

    # The published formula is cs[0] + ln(cs[1] * ln(1 - M)) / cs[2] with
    # mortality score M = 1 - exp(-exp(xb) * (exp(120 * gamma) - 1) / gamma).
    # Since -ln(1 - M) = exp(xb) * (exp(120 * gamma) - 1) / gamma exactly, it
    # reduces to the form below, which is linear in xb. Computing M first
    # loses all precision once M rounds to 1 for high-risk inputs.
    xb = constant + pheno
    phenotypic_age = (
        cs[0]
        + (np.log(-cs[1] * (np.exp(120 * gamma) - 1) / gamma) + xb) / cs[2]
    )

    return phenotypic_age
