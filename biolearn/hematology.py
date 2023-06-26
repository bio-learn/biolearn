import pandas as pd
import numpy as np


# Biomarker parameters from the paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940111/
def phenotypic_age(df):
    pheno_coefs = {
        "age": 0.0804,
        "LBDSALSI": -0.034,
        "LBDSCRSI": 0.0095,
        "glucose": 0.1953,
        "LBXCRP": 0.0954,
        "LBXLYPCT": -0.012,
        "LBXMCVSI": 0.0268,
        "LBXRDW": 0.3356,
        "LBXSAPSI": 0.00188,
        "LBXWBCSI": 0.0554,
    }

    constant = -19.9067
    gamma = 0.0077
    cs = [141.50225, -0.00553, 0.090165]
    df["LBXCRP"] = np.log(df["LBXCRP"])
    df["pheno"] = df.apply(
        lambda x: sum([x[c] * pheno_coefs[c] for c in pheno_coefs.keys()]),
        axis=1,
    )
    df["mortality_score"] = 1 - np.exp(
        (-np.exp(constant + df["pheno"]) * ((np.exp(120 * gamma) - 1) / gamma))
    )
    df["phenotypic_age"] = (
        cs[0] + np.log(cs[1] * np.log(1.00001 - df["mortality_score"])) / cs[2]
    )
    return df["phenotypic_age"]
