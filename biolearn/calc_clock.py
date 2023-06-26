import pandas as pd
import numpy as np


def trafo(x, adult_age=20):
    x = (x + 1) / (1 + adult_age)
    y = np.where(x <= 1, np.log(x), x - 1)
    return y


def anti_trafo(x, adult_age=20):
    y = np.where(
        x < 0, (1 + adult_age) * np.exp(x) - 1, (1 + adult_age) * x + adult_age
    )
    return y


def calc_alcohol_mccartney(DNAm, pheno=None):
    # Read in the Data
    Alcohol_CpGs = pd.read_csv("clock_tb/Alcohol.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = Alcohol_CpGs["CpGmarker"].isin(DNAm.columns).all()

    if CpGCheck:
        present = Alcohol_CpGs["CpGmarker"].isin(DNAm.columns)

        betas = DNAm.loc[:, Alcohol_CpGs.loc[present, "CpGmarker"]]

        tt = betas.multiply(
            Alcohol_CpGs.loc[present, "CoefficientTraining"].values, axis=1
        ).sum(axis=1)

        if pheno is None:
            return tt
        else:
            pheno["Alcohol_McCartney"] = tt
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")


def calc_bmi_mccartney(DNAm, pheno=None):
    # Read in the Data
    BMI_CpGs = pd.read_csv("clock_tb/BMI.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = BMI_CpGs["CpGmarker"].isin(DNAm.columns).all()

    if CpGCheck:
        present = BMI_CpGs["CpGmarker"].isin(DNAm.columns)

        betas = DNAm.loc[:, BMI_CpGs.loc[present, "CpGmarker"]]

        tt = betas.multiply(
            BMI_CpGs.loc[present, "CoefficientTraining"].values, axis=1
        ).sum(axis=1)

        if pheno is None:
            return tt
        else:
            pheno["BMI_McCartney"] = tt
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")


def calc_bohlin(DNAm, pheno=None):
    # Read in the Data
    Bohlin_CpGs = pd.read_csv("clock_tb/Bohlin.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = Bohlin_CpGs["CpGmarker"].isin(DNAm.columns).all()

    if CpGCheck:
        present = Bohlin_CpGs["CpGmarker"].isin(DNAm.columns)

        betas = DNAm.loc[:, Bohlin_CpGs.loc[present, "CpGmarker"]]

        tt = (
            betas.multiply(
                Bohlin_CpGs.loc[present, "Coefficient"].values, axis=1
            ).sum(axis=1)
            + 277.2421
        )

        if pheno is None:
            return tt
        else:
            pheno["PhenoAge"] = tt
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")


def calc_DNAmClockCortical(
    DNAm, pheno=None, CpGImputation=None, imputation=False
):
    # Read in the Data
    DNAmClockCortical_CpGs = pd.read_csv("clock_tb/DNAmClockCortical.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = DNAmClockCortical_CpGs["CpGmarker"].isin(DNAm.columns).all()

    if CpGCheck:
        present = DNAmClockCortical_CpGs["CpGmarker"].isin(DNAm.columns)

        betas = DNAm.loc[:, DNAmClockCortical_CpGs.loc[present, "CpGmarker"]]

        tt = (
            betas.multiply(
                DNAmClockCortical_CpGs.loc[present, "Coefficient"].values,
                axis=1,
            ).sum(axis=1)
            + 0.577682570446177
        )
        tt = tt.apply(anti_trafo)

        if pheno is None:
            return tt
        else:
            pheno["DNAmClockCortical"] = tt
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")


def calc_DNAmTL(DNAm, pheno=None, CpGImputation=None, imputation=False):
    # Read in the Data
    DNAmTL_CpGs = pd.read_csv("clock_tb/DNAmTL.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = DNAmTL_CpGs["ID"].isin(DNAm.columns).all()

    if CpGCheck:
        present = DNAmTL_CpGs["ID"].isin(DNAm.columns)

        betas = DNAm.loc[:, DNAmTL_CpGs.loc[present, "ID"]]

        tt = (
            betas.multiply(
                DNAmTL_CpGs.loc[present, "Coef"].values, axis=1
            ).sum(axis=1)
            - 7.924780053
        )

        if pheno is None:
            return tt
        else:
            pheno["DNAmTL"] = tt
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")


def calc_Hannum(DNAm, pheno=None, CpGImputation=None, imputation=False):
    # Read in the Data
    Hannum_CpGs = pd.read_csv("clock_tb/Hannum.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = Hannum_CpGs["Marker"].isin(DNAm.columns).all()

    if CpGCheck:
        present = Hannum_CpGs["Marker"].isin(DNAm.columns)

        betas = DNAm.loc[:, Hannum_CpGs.loc[present, "Marker"]]

        tt = betas.multiply(
            Hannum_CpGs.loc[present, "Coefficient"].values, axis=1
        ).sum(axis=1)

        if pheno is None:
            return tt
        else:
            pheno["Hannum"] = tt
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")


def calc_Horvath1(DNAm, pheno=None, CpGImputation=None, imputation=False):
    # Read in the Data
    Horvath1_CpGs = pd.read_csv("clock_tb/Horvath1.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = Horvath1_CpGs["CpGmarker"].isin(DNAm.columns).all()

    if CpGCheck:
        present = Horvath1_CpGs["CpGmarker"].isin(DNAm.columns)

        betas = DNAm.loc[:, Horvath1_CpGs.loc[present, "CpGmarker"]]

        tt = betas.multiply(
            Horvath1_CpGs.loc[present, "CoefficientTraining"].values, axis=1
        ).sum(axis=1)

        Horvath1 = anti_trafo(tt + 0.696)
        if pheno is None:
            return Horvath1
        else:
            pheno["Horvath1"] = Horvath1
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")


def calc_Horvath2(DNAm, pheno=None, CpGImputation=None, imputation=False):
    # Read in the Data
    Horvath2_CpGs = pd.read_csv("clock_tb/Horvath2.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = Horvath2_CpGs["CpGmarker"].isin(DNAm.columns).all()

    if CpGCheck:
        present = Horvath2_CpGs["CpGmarker"].isin(DNAm.columns)

        betas = DNAm.loc[:, Horvath2_CpGs.loc[present, "CpGmarker"]]

        tt = betas.multiply(
            Horvath2_CpGs.loc[present, "CoefficientTraining"].values, axis=1
        ).sum(axis=1)

        Horvath2 = anti_trafo(tt - 0.447119319)

        if pheno is None:
            return Horvath2
        else:
            pheno["Horvath2"] = Horvath2
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")


def calc_HRSInChPhenoAge(
    DNAm, pheno=None, CpGImputation=None, imputation=False
):
    # Read in the Data
    HRSInCHPhenoAge_CpGs = pd.read_csv("clock_tb/HRSInCHPhenoAge.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = HRSInCHPhenoAge_CpGs["CpGmarker"].isin(DNAm.columns).all()

    if CpGCheck:
        present = HRSInCHPhenoAge_CpGs["CpGmarker"].isin(DNAm.columns)

        betas = DNAm.loc[:, HRSInCHPhenoAge_CpGs.loc[present, "CpGmarker"]]

        tt = (
            betas.multiply(
                HRSInCHPhenoAge_CpGs.loc[
                    present, "CoefficientTraining"
                ].values,
                axis=1,
            ).sum(axis=1)
            + 52.8334080
        )

        if pheno is None:
            return tt
        else:
            pheno["HRSInChPhenoAge"] = tt
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")


def calc_Knight(DNAm, pheno=None, CpGImputation=None, imputation=False):
    # Read in the Data
    Knight_CpGs = pd.read_csv("clock_tb/Knight.csv")

    # Check if all necessary CpGs are in the data
    CpGCheck = Knight_CpGs["CpGmarker"].isin(DNAm.columns).all()

    if CpGCheck:
        present = Knight_CpGs["CpGmarker"].isin(DNAm.columns)

        betas = DNAm.loc[:, Knight_CpGs.loc[present, "CpGmarker"]]

        tt = (
            betas.multiply(
                Knight_CpGs.loc[present, "CoefficientTraining"].values, axis=1
            ).sum(axis=1)
            + 41.7
        )

        if pheno is None:
            return tt
        else:
            pheno["Knight"] = tt
            return pheno
    else:
        raise ValueError("Necessary CpGs are missing!")
