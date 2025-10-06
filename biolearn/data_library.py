import yaml
import pandas as pd
import numpy as np
import re
import requests
import gzip
import shutil
import os
import math
import matplotlib.pyplot as plt
import seaborn as sns


from biolearn.util import cached_download, get_data_file
from biolearn.defaults import default_cache
from biolearn.cache import NoCache
from biolearn.metadata import _iter_library_items
from io import BytesIO


def parse_after_colon(s):
    """Extract and return the substring after the first colon."""

    if isinstance(s, str) and s.strip() and ":" in s:
        return s.split(":")[1].strip()
    else:
        return s


def sex_parser(s):
    from biolearn.metadata import standardize_sex

    return standardize_sex(s)


def extract_numeric(s):
    """Extract the first numeric value from a string."""
    match = re.search(r"(\d+(\.\d+)?)", s)
    return float(match.group(1)) if match else None


def extract_informal_age(char):
    if "age (yrs)" in char:
        return extract_numeric(char["age (yrs)"])

    return None


parsers = {
    "numeric": lambda s: extract_numeric(parse_after_colon(s)),
    "string": lambda s: parse_after_colon(s),
    "sex": lambda s: sex_parser(parse_after_colon(s)),
}


def build_column_mapping(matrix_file_path, from_key_line, to_key_line):
    # Use the key line for the mapping
    mapping_df = pd.read_table(
        matrix_file_path,
        index_col=0,
        skiprows=lambda x: x != from_key_line - 1 and x != to_key_line - 1,
    )
    column_mapping = mapping_df.to_dict("records")[0]

    # Reverse the mapping if needed as key is based on first line loaded
    reverse_mapping = to_key_line < from_key_line
    if reverse_mapping:
        column_mapping = {v: k for k, v in column_mapping.items()}
    return column_mapping


def map_and_prune_columns(data, column_mapping):
    data = data.rename(columns=column_mapping)
    data = data[
        [col for col in data.columns if col in column_mapping.values()]
    ]
    return data


def load_geo_metadata(metadata_file, filekey, id_row):
    load_list = [(key, filekey[key]["row"] - 1) for key in filekey.keys()]
    load_list.sort(key=lambda x: x[1])
    load_rows = [x[1] for x in load_list]
    column_names = [x[0] for x in load_list]
    metadata = pd.read_table(
        metadata_file,
        index_col=0,
        skiprows=lambda x: x != id_row - 1 and x not in load_rows,
    )
    metadata.index = column_names
    metadata = metadata.transpose()
    metadata.index.name = "id"
    for col in metadata.columns:
        parser_name = filekey[col]["parse"]
        parser = parsers[parser_name]
        metadata[col] = metadata[col].apply(parser)
    return metadata


class QualityReport:
    def __init__(self, sample_report, methylation_site_report, summary):
        self.sample_report = sample_report
        self.methylation_site_report = methylation_site_report
        self.summary = summary

    def show(self):
        """
        Pretty prints the summary of the quality report for display, along with relevant notes based on the data.
        """
        print("Quality Report Summary\n")
        print("------------------------------------------------")

        total_samples = self.summary["total_samples"]
        total_sites = self.summary["methylation_sites"]
        missing_data_points = self.summary["missing_data_points"]

        total_data_points = total_samples * total_sites
        missing_data_percentage = (
            (missing_data_points / total_data_points) * 100
            if total_data_points > 0
            else 0
        )
        print(f"Sample Count: {total_samples}")
        print(f"Methylation Sites: {total_sites}")
        print(
            f"Missing Methylation Data: {missing_data_points} ({missing_data_percentage:.2f}%)"
        )

        high_deviation_count = self.summary["samples_with_high_deviation"]
        high_deviation_percentage = (
            high_deviation_count / total_samples
        ) * 100
        print(
            f"Samples With High Deviation: {high_deviation_count} ({high_deviation_percentage:.2f}%)"
        )

        sites_with_missing_data = self.summary[
            "sites_with_over_20_percent_missing"
        ]
        missing_data_sites_percentage = (
            sites_with_missing_data / total_sites
        ) * 100
        print(
            f"Methylation Sites With Over 20% of Reads Missing: {sites_with_missing_data} ({missing_data_sites_percentage:.2f}%)"
        )

        print("\nNotes:")
        print("------------------------------------------------")
        if missing_data_points == 0:
            print(
                "- No missing data points implies that this data has already gone through an imputation process or that low quality reads were included."
            )

        if sites_with_missing_data > 0:
            print(
                "- Your data set includes methylation sites that have over 20% of reads missing. Default imputation may replace the values for all reads from this site with a gold standard."
            )

        if high_deviation_count > 0:
            print(
                "- Your data set includes samples with a high deviation. It is likely that the methylation data for these samples has been distorted due to technical issues."
            )
        print("\n")


class GeoData:
    """
    Represents genomic data with a focus on metadata and methylation data.

    GeoData facilitates the organization and access to metadata and methylation data.

    Attributes:
        metadata (DataFrame): A pandas DataFrame where rows represent different samples
                              and columns represent different data fields.
        dnam (DataFrame): A pandas DataFrame where columns represent different samples
                          and rows represent different methylation sites.
    """

    def __init__(
        self,
        metadata,
        dnam=None,
        rna=None,
        protein_alamar=None,
        protein_olink=None,
    ):
        """
        Initializes the GeoData instance.

        Args:
            metadata (DataFrame): Metadata associated with genomic samples.
            dnam (DataFrame): Methylation data associated with genomic samples.
        """
        self.metadata = metadata
        self.dnam = dnam
        self.rna = rna
        self.protein_alamar = protein_alamar
        self.protein_olink = protein_olink

    def copy(self):
        """
        Creates a deep copy of the GeoData instance.

        Returns:
            GeoData: A new instance of GeoData with copies of the metadata and dnam DataFrames.
        """
        return GeoData(
            metadata=self.metadata.copy(deep=True),
            dnam=self.dnam.copy(deep=True) if self.dnam is not None else None,
            rna=self.rna.copy(deep=True) if self.rna is not None else None,
            protein_alamar=(
                self.protein_alamar.copy(deep=True)
                if self.protein_alamar is not None
                else None
            ),
            protein_olink=(
                self.protein_olink.copy(deep=True)
                if self.protein_olink is not None
                else None
            ),
        )

    def quality_report(self, sites=None):
        """
        Generates a quality control report for the genomic data, optionally filtered by specified methylation sites,
        and includes a detailed section reporting the missing percentage for each methylation site.

        Args:
            sites (list, optional): A list of methylation site identifiers to include in the report.
                                    If None, all sites are included.

        Returns:
            QualityReport: An object containing both detailed methylation data, a summary,
                           and a detailed section for missing percentages per site.
        """
        if self.dnam is None or self.dnam.empty:
            raise ValueError(
                "This dataset does not have Methylation data. Only methylation data is currently supported by quality_report"
            )

        methylation_data = self.dnam.copy()

        # Filter methylation data if specific sites are provided
        if sites is not None:
            methylation_data = methylation_data.loc[sites]

        # Calculate mean absolute deviation
        methylation_medians = methylation_data.median(axis=1)
        deviations = methylation_data.sub(methylation_medians, axis=0)
        mean_abs_deviation = deviations.abs().mean()
        sample_report = mean_abs_deviation.to_frame(name="deviation")

        # Include 'age' in the detailed report if relevant and if analyzing all sites
        if "age" in self.metadata.columns and sites is None:
            sample_report["age"] = self.metadata["age"]

        # Calculate missing percentage for each methylation site
        missing_percentage = methylation_data.isna().mean(axis=1)
        methylation_site_report = missing_percentage.to_frame(name="missing")

        # Summary calculations
        total_nan = methylation_data.isna().sum().sum()
        sites_with_20_percent_nan = (
            methylation_data.isna().mean(axis=1) >= 0.2
        ).sum()

        summary = {
            "missing_data_points": total_nan,
            "sites_with_over_20_percent_missing": sites_with_20_percent_nan,
            "samples_with_high_deviation": (mean_abs_deviation > 0.04).sum(),
            "total_samples": methylation_data.shape[1],
            "methylation_sites": methylation_data.shape[0],
        }

        return QualityReport(sample_report, methylation_site_report, summary)

    @classmethod
    def from_methylation_matrix(cls, matrix):
        """
        Creates a GeoData instance from a methylation matrix which can be
        either a DataFrame directly or a path to a CSV file.

        Args:
            matrix (Union[str, DataFrame]): Methylation matrix as a DataFrame or the path to the CSV file containing the matrix.

        Returns:
            GeoData: An instance of GeoData with the methylation data loaded and metadata initialized.
        """
        if isinstance(matrix, str):
            # If the input is a string, assume it's a file path and read the CSV
            dnam = pd.read_csv(
                matrix, index_col=0, na_values=["NaN", "", " "]
            ).apply(pd.to_numeric, errors="coerce")
        elif isinstance(matrix, pd.DataFrame):
            # If the input is already a DataFrame, use it directly
            dnam = matrix.copy()
        else:
            raise ValueError(
                "The matrix must be either a DataFrame or a file path to a CSV."
            )

        # Strip extra spaces from index and ensure they are strings
        dnam.index = dnam.index.str.strip().str.split("_").str[0]

        # Combine duplicate rows by averaging their values (ignoring NaN)
        dnam = dnam.groupby(dnam.index).mean()

        # Create an empty DataFrame for metadata with row identifiers as columns
        metadata = pd.DataFrame(index=dnam.columns)

        return cls(metadata, dnam)

    def save_csv(self, folder_path, name):
        """
        Saves the GeoData instance to CSV files according to the DNA Methylation Array Data Standard V-2410.

        Args:
            folder_path (str): The directory where the files will be saved.
            name (str): The base name for the saved files.

        Returns:
            None
        """
        os.makedirs(folder_path, exist_ok=True)

        # Save metadata
        if self.metadata is not None:
            metadata_to_save = self.metadata.copy()
            metadata_file = os.path.join(folder_path, f"{name}_metadata.csv")
            metadata_to_save.to_csv(metadata_file)

        # Save methylation data (dnam)
        if self.dnam is not None:
            num_samples = self.dnam.shape[1]
            if num_samples > 1000:
                parts = (num_samples - 1) // 1000 + 1
                for i in range(parts):
                    start = i * 1000
                    end = min((i + 1) * 1000, num_samples)
                    part_df = self.dnam.iloc[:, start:end]
                    file_name = os.path.join(
                        folder_path, f"{name}_methylation_part{i+1}.csv"
                    )
                    part_df.to_csv(file_name)
            else:
                file_name = os.path.join(
                    folder_path, f"{name}_methylation_part1.csv"
                )
                self.dnam.to_csv(file_name)

        # Save RNA and protein data (if present)
        if self.rna is not None:
            rna_file = os.path.join(folder_path, f"{name}_rna.csv")
            self.rna.to_csv(rna_file)
        if self.protein_alamar is not None:
            protein_file = os.path.join(
                folder_path, f"{name}_protein_alamar.csv"
            )
            self.protein_alamar.to_csv(protein_file)
        if self.protein_olink is not None:
            protein_file = os.path.join(
                folder_path, f"{name}_protein_olink.csv"
            )
            self.protein_olink.to_csv(protein_file)

    @classmethod
    def load_csv(cls, folder_path, name, series_part="all"):
        """
        Loads a GeoData instance from CSV files according to the DNA Methylation Array Data Standard V-2410.

        Args:
            folder_path (str): The directory where the files are located.
            name (str): The base name for the files.
            series_part (str or int): "all" to load all methylation parts and concatenate;
                                    otherwise, an integer specifying the part number to load.

        Returns:
            GeoData: A GeoData instance with metadata, methylation data, RNA, and protein data loaded.
        """
        # Load metadata
        metadata_file = os.path.join(folder_path, f"{name}_metadata.csv")
        if os.path.exists(metadata_file):
            metadata_df = pd.read_csv(
                metadata_file, index_col=0, keep_default_na=False
            )
        else:
            metadata_df = None

        # Load methylation data
        dnam_dfs = []
        if series_part == "all":
            files = [
                fname
                for fname in os.listdir(folder_path)
                if fname.startswith(f"{name}_methylation_part")
                and fname.endswith(".csv")
            ]
            files_sorted = sorted(
                files,
                key=lambda f: int(
                    f.split("methylation_part")[-1].split(".")[0]
                ),
            )
            for fname in files_sorted:
                part_df = pd.read_csv(
                    os.path.join(folder_path, fname),
                    index_col=0,
                    skipinitialspace=True,
                )
                dnam_dfs.append(part_df)
            dnam_df = pd.concat(dnam_dfs, axis=1) if dnam_dfs else None
        else:
            try:
                part_number = int(series_part)
            except Exception:
                raise ValueError(
                    "series_part must be 'all' or an integer value"
                )
            fname = f"{name}_methylation_part{part_number}.csv"
            file_path = os.path.join(folder_path, fname)
            dnam_df = (
                pd.read_csv(file_path, index_col=0, skipinitialspace=True)
                if os.path.exists(file_path)
                else None
            )

        # Load RNA and protein data (if available)
        rna_file = os.path.join(folder_path, f"{name}_rna.csv")
        rna_df = (
            pd.read_csv(rna_file, index_col=0, skipinitialspace=True)
            if os.path.exists(rna_file)
            else None
        )

        # In case someone saved with legacy code
        unspecified_proteomic_file = os.path.join(
            folder_path, f"{name}_protein.csv"
        )
        if os.path.exists(unspecified_proteomic_file):
            raise Exception(
                "Unspecified source proteomic file found. Please rename to specified source (e.g. protein_alamar or protein_olink) before saving."
            )
        protein_alamar_file = os.path.join(
            folder_path, f"{name}_protein_alamar.csv"
        )
        protein_alamar_df = (
            pd.read_csv(
                protein_alamar_file, index_col=0, skipinitialspace=True
            )
            if os.path.exists(protein_alamar_file)
            else None
        )

        protein_olink_file = os.path.join(
            folder_path, f"{name}_protein_olink.csv"
        )
        protein_olink_df = (
            pd.read_csv(protein_olink_file, index_col=0, skipinitialspace=True)
            if os.path.exists(protein_olink_file)
            else None
        )

        return cls(
            metadata_df,
            dnam=dnam_df,
            rna=rna_df,
            protein_alamar=protein_alamar_df,
            protein_olink=protein_olink_df,
        )


class JenAgeCustomParser:
    def __init__(self, data):
        self.data = data

    def parse(self, _):
        # Load up part 1
        jen1_meta_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103232/matrix/GSE103232_series_matrix.txt.gz"
        jen1_meta = pd.read_table(jen1_meta_url, index_col=0, skiprows=32)
        jen1_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103232/suppl/GSE103232%5Fhs%5Fblood%5Fbatch2%5Fcounts%5Frpkm.xls.gz"
        response1 = requests.get(jen1_url)
        with gzip.open(BytesIO(response1.content)) as f:
            jen1 = pd.read_excel(
                f, index_col=0, engine="xlrd"
            )  # Specify engine

        jen1 = jen1.drop(
            ["external_gene_id", "description", "gene_biotype"], axis=1
        )
        jen1 = jen1[jen1.sum(1) > 0].T

        # Age is in a 5 year range so we add 2 to the low number to get closer to the average
        age = jen1_meta.iloc[8].str[5:7].astype(int) + 2
        age.index = age.index.str[-3:]
        jen1 = jen1.join(age.rename("age"))

        # Download and read the metadata for jen2
        jen2_meta_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75337/matrix/GSE75337_series_matrix.txt.gz"
        jen2_meta = pd.read_table(jen2_meta_url, index_col=0, skiprows=34)

        # Load up part 2
        jen2_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75337/suppl/GSE75337%5Fcounts%5FRPKMs.xls.gz"
        response2 = requests.get(jen2_url)
        with gzip.open(BytesIO(response2.content)) as f:
            jen2 = pd.read_excel(
                f, index_col=0, engine="xlrd"
            )  # Specify engine

        jen2.columns = jen2.columns.str.strip(" ")
        jen2 = jen2.drop(
            ["external_gene_id", "description", "gene_biotype"], axis=1
        )
        jen2 = jen2[[c for c in jen2.columns if "blood" in c]]
        jen2 = jen2[jen2.sum(1) > 0].T

        # Age is in a 5 year range so we add 2 to the low number to get closer to the average
        age2 = jen2_meta.iloc[8].str[5:7].astype(int) + 2
        jen2 = jen2.merge(
            age2.rename("age"), left_index=True, right_index=True
        )

        # Build GeoData
        jen = pd.concat([jen1, jen2], ignore_index=True, join="inner").copy()
        metadata = jen[["age"]].copy()
        rna = jen.drop(["age"], axis=1).T.copy()
        return GeoData(metadata=metadata, rna=rna)


class ChallengeDataParser:
    def __init__(self, data):
        if data.get("id-row") is None:
            raise ValueError("Parser not valid: missing id-row")
        self.id_row = data.get("id-row")
        self.metadata = data.get("metadata")
        self.matrix_file = data.get("matrix-file")
        self.matrix_file_key_line = data.get("matrix-file-key-line")
        self.data_type = data.get("data-type")
        self.protein_matrix_url = "https://storage.googleapis.com/boa-challenge-2024/challenge_alamar_data.csv"
        self.metadata_url = "https://storage.googleapis.com/boa-challenge-2024/challenge_proteomic_metadata.csv"
        self.id_map_file = get_data_file("reference/challenge_id_map.csv")

    def parse(self, file_path):
        print("Note: This dataset will take a few minutes to load")
        # Load methylation data and metadata from GEO
        metadata = load_geo_metadata(file_path, self.metadata, self.id_row)
        dnam_data = pd.read_csv(self.matrix_file, index_col=0)
        column_mapping = build_column_mapping(
            file_path, self.matrix_file_key_line, self.id_row
        )
        fixed_dnam = map_and_prune_columns(dnam_data, column_mapping)
        geodata = GeoData.from_methylation_matrix(fixed_dnam)
        geodata.metadata = metadata

        # Load proteomic data and metadata from google cloud
        protein_matrix = pd.read_csv(self.protein_matrix_url, index_col=0)
        proteomic_metadata = pd.read_csv(self.metadata_url, index_col=0)
        proteomic_metadata["sex"] = proteomic_metadata["sex"].apply(sex_parser)

        # Use ID mapping to unify the ids
        id_map = pd.read_csv(self.id_map_file)
        plasma_to_geo_mapping = id_map.set_index("PlasmaID")["GeoID"].to_dict()
        protein_matrix.rename(columns=plasma_to_geo_mapping, inplace=True)
        proteomic_metadata.rename(index=plasma_to_geo_mapping, inplace=True)

        # Add the proteomic metadata and merge the metadata
        geodata.protein_alamar = protein_matrix.T
        merged_metadata = metadata.combine_first(proteomic_metadata)
        geodata.metadata = merged_metadata

        return geodata


class GisbyOlinkParser:
    def __init__(self, data):
        if not data.get("protein-matrix-url"):
            raise ValueError("Parser not valid: missing protein-matrix-url")
        if not data.get("metadata-url"):
            raise ValueError("Parser not valid: missing metadata-url")

        self.protein_matrix_url = data.get("protein-matrix-url")
        self.metadata_url = data.get("metadata-url")
        self.id_row = data.get("id-row", "SampleID")

    def parse(self, file_path=None):
        protein_long = pd.read_csv(self.protein_matrix_url)
        metadata = pd.read_csv(self.metadata_url, index_col=self.id_row)

        # Handle duplicate genes across panels
        gene_panel_counts = protein_long.groupby("GeneID")["Panel"].nunique()
        duplicated_genes = gene_panel_counts[
            gene_panel_counts > 1
        ].index.tolist()

        if duplicated_genes:
            protein_long["Gene_Unique"] = protein_long["GeneID"]
            mask = protein_long["GeneID"].isin(duplicated_genes)
            protein_long.loc[mask, "Gene_Unique"] = (
                protein_long.loc[mask, "GeneID"]
                + "_"
                + protein_long.loc[mask, "Panel"]
            )
            gene_column = "Gene_Unique"
        else:
            gene_column = "GeneID"

        # Create protein matrix
        protein_matrix = protein_long.pivot_table(
            index="SampleID",
            columns=gene_column,
            values="NPX",
            aggfunc="first",
        )

        # Process metadata
        if "Sex" in metadata.columns:
            metadata["sex"] = metadata["Sex"].apply(sex_parser)
            metadata.drop("Sex", axis=1, inplace=True)

        if "Age" in metadata.columns:
            metadata["age_range"] = (
                metadata["Age"].str.strip("()[]").str.replace(",", "-")
            )
            metadata.drop("Age", axis=1, inplace=True)

        for col in ["Ever_Admitted", "Fatal_Disease"]:
            if col in metadata.columns:
                metadata[col] = metadata[col].astype("bool")

        # Align samples
        common_samples = protein_matrix.index.intersection(metadata.index)
        protein_matrix = protein_matrix.loc[common_samples]
        metadata = metadata.loc[common_samples]

        return GeoData(metadata, protein_olink=protein_matrix)


class FilbinOlinkParser:
    def __init__(self, data):
        """
        Initialize the FilbinOlinkParser with paths to three Excel files.

        Args:
            data (dict): Dictionary containing:
                - metadata-path: Path to metadata Excel file
                - proteomics-path: Path to proteomics Excel file
                - mappings-path: Path to olink mappings Excel file
        """
        required_paths = ["metadata-path", "proteomics-path", "mappings-path"]
        for path_key in required_paths:
            if not data.get(path_key):
                raise ValueError(f"Parser not valid: missing {path_key}")

        self.metadata_path = data["metadata-path"]
        self.proteomics_path = data["proteomics-path"]
        self.mappings_path = data["mappings-path"]

        self.age_mapping = {
            1: "20-34",
            2: "36-49",
            3: "50-64",
            4: "65-79",
            5: "80+",
        }

    def parse(self, file_path=None):
        """
        Parse the three Excel files and create a GeoData object.

        Returns:
            GeoData: Processed data in GeoData format
        """
        metadata_df = (
            pd.read_excel(
                self.metadata_path,
                sheet_name=0,
                usecols=["Public ID", "Age cat"],
            )
            .rename(columns={"Public ID": "PublicID", "Age cat": "AgeCat"})
            .set_index("PublicID")
        )

        proteomics_df = pd.read_excel(self.proteomics_path)
        proteomics_df = proteomics_df.drop(columns=["Day"], errors="ignore")

        id_col = (
            "Public ID"
            if "Public ID" in proteomics_df.columns
            else proteomics_df.columns[0]
        )
        proteomics_df = proteomics_df.rename(
            columns={id_col: "SampleID"}
        ).set_index("SampleID")

        mappings_df = pd.read_excel(
            self.mappings_path,
            sheet_name="2A-Olink-Assay",
            usecols=["OlinkID", "Assay"],
            skiprows=1,
        )

        # Create mapping and track gene-to-OlinkIDs for duplicates
        olink_to_gene = dict(zip(mappings_df["OlinkID"], mappings_df["Assay"]))
        gene_to_olinks = {}
        for olink_id, gene in zip(
            mappings_df["OlinkID"], mappings_df["Assay"]
        ):
            gene_to_olinks.setdefault(gene, []).append(olink_id)

        # Rename columns and track which OlinkIDs are in the data
        column_mapping = {
            col: olink_to_gene.get(col, col) for col in proteomics_df.columns
        }
        proteomics_df = proteomics_df.rename(columns=column_mapping)

        # Handle duplicate columns after mapping
        duplicate_cols = proteomics_df.columns[
            proteomics_df.columns.duplicated()
        ].unique()
        if len(duplicate_cols) > 0:
            # TODO: Figure out why olink provides duplicate IDs for the same protein
            proteomics_df = proteomics_df.loc[
                :, ~proteomics_df.columns.duplicated(keep="first")
            ]

        # Expand metadata for each sample
        expanded_metadata_rows = []
        for sample_id in proteomics_df.index:
            sample_str = str(sample_id)

            # Extract metadata ID (numeric part before underscore)
            try:
                metadata_id = (
                    int(sample_str.split("_")[0])
                    if "_" in sample_str
                    else int(sample_str)
                )
            except (ValueError, IndexError):
                continue

            if metadata_id in metadata_df.index:
                meta_row = metadata_df.loc[metadata_id].copy()
                meta_row.name = sample_id
                expanded_metadata_rows.append(meta_row)

        expanded_metadata = (
            pd.DataFrame(expanded_metadata_rows)
            if expanded_metadata_rows
            else metadata_df.copy()
        )

        # Map age categories to age ranges
        if "AgeCat" in expanded_metadata.columns:
            expanded_metadata["age_range"] = expanded_metadata["AgeCat"].map(
                self.age_mapping
            )
            expanded_metadata = expanded_metadata.drop(columns=["AgeCat"])

        # Align samples
        common_samples = proteomics_df.index.intersection(
            expanded_metadata.index
        )
        if len(common_samples) == 0:
            print(
                "Warning: No common samples found. Using all proteomics samples."
            )
            common_samples = proteomics_df.index

        proteomics_aligned = proteomics_df.loc[common_samples]
        metadata_aligned = expanded_metadata.loc[common_samples]

        return GeoData(metadata_aligned, protein_olink=proteomics_aligned)


class GeoMatrixParser:
    seperators = {"space": " ", "comma": ",", "tab": "\t"}

    def __init__(self, data):
        if data.get("id-row") is None:
            raise ValueError("Parser not valid: missing id-row")
        self.id_row = data.get("id-row")
        self.metadata = data.get("metadata")
        self.matrix_start = data.get("matrix-start")
        self.matrix_file = data.get("matrix-file")
        self.matrix_file_seperator = self.seperators.get(
            data.get("matrix-file-seperator")
        )
        self.matrix_file_key_line = data.get("matrix-file-key-line")
        self.matrix_file_format = data.get("matrix-file-format")
        self.data_type = data.get("data-type")

    def parse(self, file_path):
        metadata = load_geo_metadata(file_path, self.metadata, self.id_row)
        if self.matrix_start:
            matrix_data = pd.read_table(
                file_path, index_col=0, skiprows=self.matrix_start - 1
            )
            matrix_data = matrix_data.drop(
                ["!series_matrix_table_end"], axis=0
            )
            matrix_data.index.name = "id"
        elif self.matrix_file:
            matrix_file_path = cached_download(self.matrix_file)
            print(
                "Note: This dataset will take a few minutes to load and may fail if you have insufficient memory"
            )
            df = pd.read_csv(
                matrix_file_path, index_col=0, sep=self.matrix_file_seperator
            )
            if self.matrix_file_format == "pvalue":
                reading_df = df.iloc[:, ::2]
                pval_df = df.iloc[:, 1::2]
                pval_df = pval_df.replace("<1E-16", "0", regex=False).astype(
                    float, errors="ignore"
                )
                pval_df = pval_df.map(lambda x: np.nan if x > 0.05 else 0)
                # NaN values in pval_df will cause corresponding values in methylation_df to be NaN
                matrix_data = reading_df + pval_df.values
                matrix_data = self._remap_and_prune_columns(
                    matrix_data, file_path
                )

            elif self.matrix_file_format == "standard":
                matrix_data = df
                matrix_data = self._remap_and_prune_columns(
                    matrix_data, file_path
                )

            else:
                raise ValueError(
                    f"Unsupported matrix file format: {self.matrix_file_format}"
                )

        if self.data_type == "rna":
            return GeoData(metadata, rna=matrix_data)
        else:
            return GeoData(metadata, dnam=matrix_data)

    def _remap_and_prune_columns(self, data, matrix_file_path):
        if self.matrix_file_key_line is None:
            # No key line for mapping so assume the ordering is sufficient
            header_row = pd.read_table(
                matrix_file_path,
                index_col=0,
                header=None,
                skiprows=lambda x: x != self.id_row - 1,
                nrows=1,
            )
            column_mapping = dict(zip(data.columns, header_row.iloc[0]))
        else:
            # Use the key line for the mapping
            mapping_df = pd.read_table(
                matrix_file_path,
                index_col=0,
                skiprows=lambda x: x != self.id_row - 1
                and x != self.matrix_file_key_line - 1,
            )
            column_mapping = mapping_df.to_dict("records")[0]

            # Reverse the mapping if needed as key is based on first line loaded
            reverse_mapping = self.id_row < self.matrix_file_key_line
            if reverse_mapping:
                column_mapping = {v: k for k, v in column_mapping.items()}

        data = data.rename(columns=column_mapping)
        data = data[
            [col for col in data.columns if col in column_mapping.values()]
        ]
        return data

    def _metadata_load_list(self):
        load_list = [
            (key, self.metadata[key]["row"] - 1)
            for key in self.metadata.keys()
        ]
        load_list.sort(key=lambda x: x[1])
        return load_list


class AutoScanGeoMatrixParser:
    """
    Retrieve Series Metadata via geo2r API

    When dealing with a large number of geo series data, instead of downloading the matrix file for each series to read metadata information, we can reduce the time spent generating the library.yaml by using the geo2r API to obtain metadata during data loading.
    """

    def __init__(self, data) -> None:
        self._check_keys_exist(
            data, ["matrix_file", "metadata_keys_parse", "metadata_query"]
        )
        self.matrix_file = data["matrix_file"]
        self.metadata_keys_parse = data["metadata_keys_parse"]
        self.metadata_keys = list(self.metadata_keys_parse.keys())
        self.metadata_query_url = data["metadata_query"]

    def _check_keys_exist(self, data, keys: list[str]):
        for key in keys:
            if key not in data:
                raise ValueError(f"Parser not valid: missing {key}")

    def _create_metadata(self, metadata_query_url):
        response = requests.get(metadata_query_url)
        json_data = response.json()
        smaples = [item for item in json_data["GeoMetaData"]]
        return self._convert_to_metadata_df(self.metadata_keys, smaples)

    def _map_characteristics_to_dict(self, sample):

        def extract_key(item):
            if "tag" in item:
                return item["tag"]
            else:
                return item["value"].split(":")[0]

        def extract_value(item):
            if "tag" in item:
                return item["value"]
            else:
                return item["value"].split(":")[1]

        characteristics = sample["entity"]["sample"]["channels"][0][
            "characteristics"
        ]

        return {
            extract_key(char).lower(): extract_value(char)
            for char in characteristics
        }

    def _convert_characteristics_to_df_cols_data(self, metadata_keys, item):
        cols_data = []

        characteristics = self._map_characteristics_to_dict(item)

        for key in metadata_keys:
            parse_type = self.metadata_keys_parse[key]
            if parse_type == "sex":
                cols_data.append(sex_parser(characteristics.get(key)))
            elif parse_type == "numeric":

                if key == "age":
                    if key in characteristics:
                        cols_data.append(extract_numeric(characteristics[key]))
                    else:
                        cols_data.append(extract_informal_age(characteristics))
                else:
                    cols_data.append(
                        extract_numeric(characteristics.get(key, ""))
                    )
            else:
                cols_data.append(characteristics.get(key))

        return cols_data

    def _convert_to_metadata_df(self, metadata_keys, samples):

        data = {
            sample["acc"]: self._convert_characteristics_to_df_cols_data(
                metadata_keys, sample
            )
            for sample in samples
        }

        df = pd.DataFrame(data).transpose()
        df = df.reset_index(drop=False)
        columns = ["id"]
        columns.extend(metadata_keys)
        df.columns = columns
        return df

    def _get_matrix_file_size_from_url(self, url):
        try:
            response = requests.head(url)
            response.raise_for_status()

            content_length = response.headers.get("Content-Length")

            if content_length is None:
                print("Content-Length header is not present.")
                return None

            file_size = int(content_length)
            return file_size

        except requests.RequestException as e:
            print(f"HTTP request failed: {e}")
            return None

    def _create_matrix(self, matrix_file):
        matrix_file_path = cached_download(matrix_file)

        row_num = self._get_matrix_table_row_num(matrix_file_path)

        matrix_data = pd.read_table(
            matrix_file_path, index_col=0, skiprows=row_num
        )
        matrix_data = matrix_data.drop(["!series_matrix_table_end"], axis=0)
        matrix_data.index.name = "id"

        return matrix_data

    def _ungzip_file(self, matrix_file, output_file):
        with gzip.open(matrix_file, "rb") as f_in:
            with open(output_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    def _find_matrix_table_begin_line_num(self, matrix_file):
        with open(matrix_file, "r", encoding="utf-8") as file:
            for line_number, line in enumerate(file, start=1):
                if "!series_matrix_table_begin" in line:
                    return line_number
        return -1

    def _get_matrix_table_row_num(self, matrix_file):

        output_file = f"{matrix_file}.txt"
        try:
            self._ungzip_file(matrix_file, output_file)
            line_num = self._find_matrix_table_begin_line_num(output_file)
        finally:
            if os.path.exists(output_file):
                os.remove(output_file)

        return line_num

    def parse(self, _):
        """
        Parse metadata and matrix data for a geo series.

        This function retrieves metadata and matrix data, raising an error if the matrix is empty or None.
        """
        metadata = self._create_metadata(self.metadata_query_url)
        matrix = self._create_matrix(self.matrix_file)
        if matrix is None or matrix.empty:
            raise NoMatrixDataError(
                f"Series {metadata['id']} has no matrix data"
            )

        return GeoData(metadata, matrix)


class NoMatrixDataError(Exception):

    def __init__(self, message):
        super().__init__(message)


class DataSource:
    """
    Represents a single data source in the DataLibrary.

    This class encapsulates the details of a data source including metadata about the source and
    functionality to load the data.

    Raises:
        ValueError: If any of the required fields ('id', 'path', 'parser') are missing
                    during initialization.
    """

    REQUIRED_FIELDS = {
        "id": "'id' key is missing in item",
        "path": "'path' key is missing in item",
        "parser": "'parser' key is missing in item",
    }

    CACHE_CATEGORY = "data_source"
    CACHE_VERSION = "v2"

    def __init__(self, source_definition, cache=None):
        """
        Initializes the DataSource instance with configuration data and an optional cache mechanism.
        This method parses a dictionary typically loaded from a YAML configuration file for a data source.
        It checks for essential fields, sets up attributes, and configures a parser for data handling.

        Args:
            source_definition (dict): A dictionary containing the data source's properties. Must include keys like 'id',
                        'path', 'parser', and optionally 'title', 'summary', 'format', and 'organism'.
            cache (object, optional): An object that adheres to the caching interface used in the caching module.
                        If no cache is provided, a default NoCache instance is used.

        Raises:
            ValueError: If any of the required fields ('id', 'path', 'parser') are missing in the input data.
        """
        for field, error_message in self.REQUIRED_FIELDS.items():
            setattr(self, field, source_definition.get(field))
            if getattr(self, field) is None:
                raise ValueError(error_message)
        self.cache = cache if cache else NoCache()
        self.source_definition = source_definition
        self.title = source_definition.get(
            "title", ""
        )  # Default empty string if title is not provided
        self.summary = source_definition.get(
            "summary", ""
        )  # Default empty string if summary is not provided
        self.format = source_definition.get(
            "format", ""
        )  # Default empty string if format is not provided
        self.organism = source_definition.get(
            "organism", ""
        )  # Default empty string if organism is not provided
        self.tags = source_definition.get(
            "tags", []
        )  # Default empty list if tags are not provided

        self.parser = self._create_parser(source_definition["parser"])

    def load(self):
        """
        Loads the data from the source.
        Returns:
            GeoData: An instance of the GeoData class containing the parsed geographical data.
        """

        if self.tags and "work_needed" in self.tags:
            self._show_work_needed_warning()

        cached = self.cache.get(
            self.id, self.CACHE_CATEGORY, self.CACHE_VERSION
        )
        if cached is not None:
            return cached

        data = self.parser.parse(self.path)
        self.cache.store(
            self.id, data, self.CACHE_CATEGORY, self.CACHE_VERSION
        )
        return data

    def __repr__(self):
        """
        Represents the DataSource instance as a string.

        Returns:
            str: A string representation of the DataSource instance.
        """
        return f"DataSource(ID: {self.id}, Title: {self.title})"

    def _create_parser(self, parser_data):
        parser_type = parser_data.get("type")
        if parser_type is None:
            raise ValueError("Parser type is missing")
        if parser_type == "jen-age-custom":
            return JenAgeCustomParser(parser_data)
        if parser_type == "geo-matrix":
            return GeoMatrixParser(parser_data)
        if parser_type == "geo-matrix-auto-scan":
            return AutoScanGeoMatrixParser(parser_data)
        if parser_type == "biomarkers-challenge-2024":
            return ChallengeDataParser(parser_data)
        if parser_type == "gisby-olink":
            return GisbyOlinkParser(parser_data)
        if parser_type == "filbin-olink":
            return FilbinOlinkParser(parser_data)
        raise ValueError(f"Unknown parser type: {parser_type}")

    def _show_work_needed_warning(self):
        """Display warning for datasets marked as needing work"""
        warning_message = (
            f"\nWARNING: Dataset {self.id} is marked as 'work_needed'.\n"
            "This indicates the dataset may have one or more issues:\n"
            "- DNA methylation data may not be populated\n"
            "- Sex information is missing or not standardized\n"
            "- Age values are not numeric or not standardized\n"
            "Please verify and standardize these fields before using in production.\n"
        )
        print(warning_message)


def parse_library_file(library_file, cache=None):
    data = yaml.safe_load(library_file)
    if data is None:
        raise ValueError("File must be YAML format with 'items' at root")
    if "items" in data:
        data_sources = [
            DataSource(item, cache if cache else NoCache())
            for item in data["items"]
        ]
        return data_sources
    else:
        raise ValueError("File must be YAML format with 'items' at root")


class DataLibrary:
    """
    Manages a collection of data sources for biomarkers research.

    The DataLibrary class is responsible for loading, storing, and retrieving data sources.
    Data sources are defined in a library file and new sources can easily be added at runtime.
    Currently DNA methylation data from GEO is supported.
    """

    def __init__(self, library_file=None, cache=None):
        """
        Initializes the DataLibrary instance with an optional library file and cache mechanism.

        Args:
            library_file (str, optional): The path to the library file. If None, the default biolearn
                                          library file is loaded.
            cache (object, optional): An object that adheres to the caching interface used in the caching module. If None, the
                                      default cache is used. This cache will be used by all returned
                                      data sources
        """
        self.cache = cache if cache else default_cache()
        self.sources = []
        if library_file is None:
            curated_library_file = get_data_file("library.yaml")
            self.load_sources(curated_library_file)
            autoscan_library = get_data_file("geo_autoscan_library.yaml")
            self.load_sources(autoscan_library)
        else:
            self.load_sources(library_file)

    def load_sources(self, library_file):
        """
        Loads data sources from a given library file appending them to the current set of data sources.

        Args:
            library_file (str): The file path of the library file to load data sources from.
        """
        with open(library_file, "r") as f:
            data_sources = self._parse_library_file(f)
            self.sources.extend(data_sources)

    def get(self, source_id):
        """
        Retrieves a data source by its identifier.

        Args:
            source_id (str): The identifier of the data source to retrieve.

        Returns:
            The data source with the given identifier if found, otherwise None.
        """
        for source in self.sources:
            if source.id == source_id:
                return source
        return None

    def lookup_sources(self, organism=None, format=None):
        """
        Looks up data sources based on the specified organism and/or format.

        Args:
            organism (str, optional): The organism to filter the data sources by.
            format (str, optional): The format to filter the data sources by.

        Returns:
            A list of data sources that match the specified organism and format criteria.
        """
        matches = []
        for source in self.sources:
            if (organism is None or source.organism == organism) and (
                format is None or source.format == format
            ):
                matches.append(source)
        return matches

    def _parse_library_file(self, library_file):
        return parse_library_file(library_file, self.cache)

    def search(self, **criteria):
        """
        Search and preview metadata across all available datasets without loading them.

        This method allows you to explore what datasets are available and their metadata
        characteristics before deciding which ones to load. It's particularly useful for
        discovering datasets that match specific criteria like sex, age, or other metadata fields.

        Parameters
        ----------
        criteria : keyword arguments
            Keyword arguments for filtering datasets.
            Common filters include:

            - sex (str): Filter by sex ("male", "female", "unknown")
            - min_age (float): Minimum age threshold
            - max_age (float): Maximum age threshold

        Returns
        -------
        pandas.DataFrame
            A DataFrame with columns including 'series_id' and available
            metadata fields for each matching dataset.

        Examples
        --------
        >>> # Find all datasets with female subjects
        >>> library = DataLibrary()
        >>> female_datasets = library.search(sex="female")

        >>> # Find datasets with elderly subjects (70+ years)
        >>> elderly_datasets = library.search(min_age=70)

        >>> # Find male datasets with subjects over 50
        >>> male_elderly = library.search(sex="male", min_age=50)

        >>> # View available metadata fields
        >>> all_datasets = library.search()
        >>> print(all_datasets.columns.tolist())

        Notes
        -----
        Sex encoding follows the DNA Methylation Array Data Standard:
        - 0 = female
        - 1 = male
        - NaN = unknown/missing
        """

        hits = []

        try:
            for sid, entry in _iter_library_items():
                # Resolve metadata dictionary - always use a copy to avoid modifying original data
                meta = {}

                # First check for direct metadata
                if "metadata" in entry:
                    meta = entry["metadata"].copy()

                # Then check for direct top-level fields (but don't overwrite existing)
                if "sex" in entry and "sex" not in meta:
                    meta["sex"] = entry["sex"]
                if "age" in entry and "age" not in meta:
                    age_val = entry["age"]
                    if isinstance(age_val, (int, float, str)):
                        try:
                            meta["age"] = float(age_val)
                        except (ValueError, TypeError):
                            pass

                # Apply sex filter
                if "sex" in criteria:
                    wanted = criteria["sex"].lower()
                    raw_sex = meta.get("sex")

                    # Check if dataset has sex information available
                    parser = entry.get("parser", {})
                    # Handle both 'metadata' and 'metadata_keys_parse' structures
                    parser_metadata = parser.get("metadata", {})
                    if not parser_metadata and "metadata_keys_parse" in parser:
                        parser_metadata = parser.get("metadata_keys_parse", {})
                    has_sex_parser = "sex" in parser_metadata

                    # If we have direct sex metadata, try to match it
                    if raw_sex is not None:
                        sex_str = None

                        # Convert to standardized string for comparison
                        if isinstance(raw_sex, str):
                            sex_lower = raw_sex.strip().lower()
                            # Map common variations
                            if sex_lower in ["f", "female"]:
                                sex_str = "female"
                            elif sex_lower in ["m", "male"]:
                                sex_str = "male"
                            elif sex_lower in ["unknown", "u", ""]:
                                sex_str = "unknown"
                            else:
                                sex_str = sex_lower
                        elif isinstance(
                            raw_sex, (int, float)
                        ) and not math.isnan(raw_sex):
                            # Handle numeric encoding
                            sex_map = {
                                0: "female",  # Standard
                                1: "male",  # Standard
                                2: "female",  # GEO encoding (reversed)
                            }
                            sex_str = sex_map.get(int(raw_sex), "unknown")
                        else:
                            sex_str = "unknown"

                        # Skip if sex doesn't match
                        if sex_str != wanted:
                            continue

                    elif not has_sex_parser:
                        # No sex information at all - skip this dataset
                        continue
                    # If has_sex_parser but no direct metadata, include the dataset
                    # because we can't determine the actual sex without loading data

                # Apply age filters
                if "min_age" in criteria:
                    age = meta.get("age")
                    if age is None or age < criteria["min_age"]:
                        continue

                if "max_age" in criteria:
                    age = meta.get("age")
                    if age is None or age > criteria["max_age"]:
                        continue

                # Create result entry with series_id
                result_entry = {"series_id": sid}

                # Add metadata fields that exist
                for key, value in meta.items():
                    if value is not None:
                        result_entry[key] = value

                hits.append(result_entry)
        except Exception:
            # If iteration fails, return empty DataFrame
            return pd.DataFrame(columns=["series_id"])

        # Always return a DataFrame, even if empty
        if not hits:
            return pd.DataFrame(columns=["series_id"])

        return pd.DataFrame(hits)
