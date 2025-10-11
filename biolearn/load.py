import pandas as pd
from biolearn.util import cached_download
from pathlib import Path
import gzip
import warnings


def verify_metadata_alignment(omics_dict, metadata):
    """
    Verify that metadata matches omics sample IDs.
    - Adds missing metadata entries.
    - Warns if metadata contains unmatched IDs.
    
    Parameters
    ----------
    omics_dict : dict[str, pd.DataFrame]
        Dictionary of omics matrices keyed by modality (e.g. "rna", "methylation").
        Each DataFrame should have sample IDs as index.
    metadata : pd.DataFrame
        Metadata with sample IDs as index.
    
    Returns
    -------
    metadata : pd.DataFrame
        Updated metadata aligned with all omics matrices.
    """

    # Collect all sample IDs from omics
    all_sample_ids = set()
    for _, df in omics_dict.items():
        all_sample_ids.update(df.index)

    # Track missing metadata
    missing_ids = all_sample_ids - set(metadata.index)
    if missing_ids:
        # Convert set to list (or sorted list) to avoid Pandas error
        new_rows = pd.DataFrame(index=list(missing_ids), columns=metadata.columns)
        metadata = pd.concat([metadata, new_rows])
        warnings.warn(
            f"{len(missing_ids)} metadata entries were missing. "
            f"Blank rows added for: {list(missing_ids)[:5]}{'...' if len(missing_ids) > 5 else ''}"
        )

    # Track extra metadata
    extra_ids = set(metadata.index) - all_sample_ids
    if extra_ids:
        warnings.warn(
            f"{len(extra_ids)} metadata entries exist without matching omics samples. "
            f"Examples: {list(extra_ids)[:5]}{'...' if len(extra_ids) > 5 else ''}"
        )

    return metadata



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


# ---------------------------------------------------------------------------
# SAS‑XPT helper that can read plain or gzipped files
# ---------------------------------------------------------------------------


def _load_sas(path, columns, name, idx="SEQN") -> pd.DataFrame:
    """Read an XPORT file (optionally gzipped) into a DataFrame and verify schema."""

    def _read(src):
        # When *src* is a file‑like object pandas needs the explicit format
        return pd.read_sas(src, format="xport", index=idx)

    path = Path(path)

    # Detect gzip via first two bytes 0x1f 0x8b
    with path.open("rb") as fh:
        magic = fh.read(2)

    if magic == b"\x1f\x8b":
        with gzip.open(path, "rb") as gz:
            df = _read(gz)
    else:
        df = _read(str(path))

    df.index = df.index.astype(int)

    if columns is not None:
        missing = [c for c in columns if c not in df.columns]
        if missing:
            raise RuntimeError(f"[{name}] missing expected cols: {missing}")
        df = df[columns]

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
    cycle = year - 1
    dem = _load_sas(
        cached_download(
            f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle}/DataFiles/DEMO_{suffix}.xpt"
        ),
        ["RIAGENDR", "RIDAGEYR"],
        name="DEMO",
    )
    gluc = _load_sas(
        cached_download(
            f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle}/DataFiles/GLU_{suffix}.xpt"
        ),
        ["LBDGLUSI"],
        name="GLU",
    )
    cbc = _load_sas(
        cached_download(
            f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle}/DataFiles/CBC_{suffix}.xpt"
        ),
        cbc_sub,
        name="CBC",
    )
    hdl = _load_sas(
        cached_download(
            f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle}/DataFiles/HDL_{suffix}.xpt"
        ),
        ["LBDHDDSI"],
        name="HDL",
    )
    bioc = _load_sas(
        cached_download(
            f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle}/DataFiles/BIOPRO_{suffix}.xpt"
        ),
        ["LBDSALSI", "LBDSCRSI", "LBXSAPSI"],
        name="BIOPRO",
    )

    if year == 2010:  # only 2010 cycle exposes CRP here
        crp = _load_sas(
            cached_download(
                f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle}/DataFiles/CRP_{suffix}.xpt"
            ),
            ["LBXCRP"],
            name="CRP",
        )

    mortality_file = cached_download(
        f"https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_{cycle}_{year}_MORT_2019_PUBLIC.dat",
    )
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
    # Wrap omics tables into dictionary for verification
    omics = {"dem": dem, "gluc": gluc, "cbc": cbc, "bioc": bioc, "hdl": hdl}
    if year == 2010:
        omics["crp"] = crp

    # Verify metadata alignment
    df = verify_metadata_alignment(omics, df)
    # rename columns
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
