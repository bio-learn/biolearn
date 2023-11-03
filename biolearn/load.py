import pandas as pd
from biolearn.util import cached_download

MG_PER_DL_TO_MMOL_PER_L = 0.05551


def load_fhs():
    """
    Loads data from the Framingham Heart Study

    Returns
    -------
    df: Pandas.Dataframe
        A pandas dataframe where each row represents an individual and each column represents a measurement about that individual
    """
    public_link = "https://raw.githubusercontent.com/singator/bdah/master/data/frmgham2.csv"
    df = pd.read_csv(
        public_link,
        index_col=0,
        usecols=[
            "RANDID",
            "PERIOD",
            "AGE",
            "SEX",
            "DEATH",
            "TIMEDTH",
            "GLUCOSE",
        ],
    )
    df = df[df["PERIOD"] == 1].drop("PERIOD", axis=1)
    df["TIMEDTH"] = df["TIMEDTH"] / 30.437  # days to months
    df.index.name = "id"

    df = df.rename(
        {
            "AGE": "age",
            "SEX": "sex",
            "GLUCOSE": "glucose",
            "DEATH": "is_dead",
            "TIMEDTH": "months_until_death",
        },
        axis=1,
    )

    # standardize glucose units
    df["glucose"] = df["glucose"].apply(lambda g: g * MG_PER_DL_TO_MMOL_PER_L)

    return df


def load_nhanes(year):
    """Loads data from the National Health and Nutrition Examination Survey

    Parameters
    ----------
    year : number
        A year number for which to load data. Not that NHANES data comes in two year groupings and the year passed in should be the later year.
        Supported inputs are 2010 and 2012

    Returns
    -------
    df: Pandas.Dataframe
        A pandas dataframe where each row represents an individual and each column represents a measurement about that individual
    """
    cbc_sub = [
        "LBXRDW",
        "LBXWBCSI",
        "LBXLYPCT",
        "LBXMCVSI",
        "LBDLYMNO",
        "LBXRBCSI",
        "LBXHGB",
        "LBXPLTSI",
        "LBXMCHSI",
        "LBXBAPCT",
    ]
    known_nhanes_year_suffix = {2010: "F", 2012: "G"}
    if not year in known_nhanes_year_suffix:
        raise ValueError(
            f"Unknown year {year}. Can only load for known available years {known_nhanes_year_suffix.keys}"
        )
    suffix = known_nhanes_year_suffix[year]
    dem_file = cached_download(
        f"https://wwwn.cdc.gov/Nchs/Nhanes/{year-1}-{year}/DEMO_{suffix}.XPT"
    )
    gluc_file = cached_download(
        f"https://wwwn.cdc.gov/Nchs/Nhanes/{year-1}-{year}/GLU_{suffix}.XPT"
    )
    cbc_file = cached_download(
        f"https://wwwn.cdc.gov/Nchs/Nhanes/{year-1}-{year}/CBC_{suffix}.XPT"
    )
    bioc_file = cached_download(
        f"https://wwwn.cdc.gov/Nchs/Nhanes/{year-1}-{year}/BIOPRO_{suffix}.XPT"
    )
    mortality_file = cached_download(
        f"https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_{year-1}_{year}_MORT_2019_PUBLIC.dat"
    )
    crp_file = cached_download(
        f"https://wwwn.cdc.gov/Nchs/Nhanes/{year-1}-{year}/CRP_{suffix}.XPT"
    )
    hdl_file = cached_download(
        f"https://wwwn.cdc.gov/Nchs/Nhanes/{year-1}-{year}/HDL_{suffix}.XPT"
    )
    dem = pd.read_sas(dem_file, index="SEQN")[["RIAGENDR", "RIDAGEYR"]]
    dem.index = dem.index.astype(int)
    gluc = pd.read_sas(gluc_file, index="SEQN")["LBDGLUSI"]
    gluc.index = gluc.index.astype(int)
    cbc = pd.read_sas(cbc_file, index="SEQN")[cbc_sub]
    cbc.index = cbc.index.astype(int)
    hdl = pd.read_sas(hdl_file, index="SEQN")["LBDHDDSI"]
    hdl.index = hdl.index.astype(int)
    # clumsy hack since 2012 doesn't have the CRP data. Will remove pending refactor of loading code
    if year == 2010:
        crp = pd.read_sas(crp_file, index="SEQN")["LBXCRP"]
        crp.index = crp.index.astype(int)
    bioc = pd.read_sas(bioc_file, index="SEQN")[
        ["LBDSALSI", "LBDSCRSI", "LBXSAPSI"]
    ]
    bioc.index = bioc.index.astype(int)
    mort = pd.read_fwf(
        mortality_file,
        index_col=0,
        header=None,
        widths=[14, 1, 1, 3, 1, 1, 1, 4, 8, 8, 3, 3],
    )
    mort.index = mort.index.rename("SEQN")
    dead = mort[mort[1] == 1][[2, 10]].astype(int)
    dead.columns = ["MORTSTAT", "PERMTH_EXM"]
    # clumsy hack since 2012 doesn't have the CRP data. Will remove pending refactor of loading code
    if year == 2010:
        df = pd.concat([dem, gluc, cbc, crp, bioc, hdl, dead], axis=1).dropna()
    else:
        df = pd.concat([dem, gluc, cbc, bioc, hdl, dead], axis=1).dropna()
    df.index.name = "id"
    df = df.rename(
        {
            "RIDAGEYR": "age",
            "RIAGENDR": "sex",
            "LBDGLUSI": "glucose",
            "MORTSTAT": "is_dead",
            "PERMTH_EXM": "months_until_death",
            "LBXWBCSI": "white_blood_cell_count",
            "LBXLYPCT": "lymphocyte_percent",
            "LBXRDW": "red_blood_cell_distribution_width",
            "LBXMCVSI": "mean_cell_volume",
            "LBDLYMNO": "lymphocyte_number",
            "LBXRBCSI": "red_blood_cell_count",
            "LBXHGB": "hemoglobin",
            "LBXPLTSI": "platelet_count",
            "LBXMCHSI": "mean_cell_hemoglobin",
            "LBXBAPCT": "basophil_percent",
            "LBDHDDSI": "hdl_cholesterol",
            "LBXCRP": "c_reactive_protein",
            "LBDSALSI": "albumin",
            "LBDSCRSI": "creatinine",
            "LBXSAPSI": "alkaline_phosphate",
        },
        axis=1,
    )
    df = df.rename({"LB2RDW": "LBXRDW", "LB2WBCSI": "LBXWBCSI"}, axis=1)
    return df
