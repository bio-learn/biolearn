import yaml
import pandas as pd
import numpy as np
import re
import requests
import gzip
import shutil
import os


from biolearn.util import cached_download, get_data_file
from biolearn.defaults import default_cache
from biolearn.cache import NoCache
from io import BytesIO


def parse_after_colon(s):
    """Extract and return the substring after the first colon."""

    if isinstance(s, str) and s.strip() and ":" in s:
        return s.split(":")[1].strip()
    else:
        return s


def sex_parser(s):
    if isinstance(s, str):
        s_lower = s.lower().strip()
        if s_lower in ["female", "f"]:
            return 1
        elif s_lower in ["male", "m"]:
            return 2
    return 0


def extract_numeric(s):
    """Extract the first numeric value from a string."""
    match = re.search(r"(\d+(\.\d+)?)", s)
    return float(match.group(1)) if match else None


def extract_informal_age(char) -> int:
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

    def __init__(self, metadata, dnam=None, rna=None, protein=None):
        """
        Initializes the GeoData instance.

        Args:
            metadata (DataFrame): Metadata associated with genomic samples.
            dnam (DataFrame): Methylation data associated with genomic samples.
        """
        self.metadata = metadata
        self.dnam = dnam
        self.rna = rna
        self.protein = protein

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
            dnam = pd.read_csv(matrix, index_col=0)
        elif isinstance(matrix, pd.DataFrame):
            # If the input is already a DataFrame, use it directly
            dnam = matrix
        else:
            raise ValueError(
                "The matrix must be either a DataFrame or a file path to a CSV."
            )

        # Process row identifiers
        dnam.index = dnam.index.str.split("_").str[0]

        # Combine duplicate rows by averaging their values
        dnam = dnam.groupby(dnam.index).mean()

        # Create an empty DataFrame for metadata with row identifiers as columns
        metadata = pd.DataFrame(index=dnam.columns)

        return cls(metadata, dnam)


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
        geodata.protein = protein_matrix
        merged_metadata = metadata.combine_first(proteomic_metadata)
        geodata.metadata = merged_metadata

        return geodata


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

        self.parser = self._create_parser(source_definition["parser"])

    def load(self):
        """
        Loads the data from the source.
        Returns:
            GeoData: An instance of the GeoData class containing the parsed geographical data.
        """
        cached = self.cache.get(self.id)
        if cached:
            return cached
        else:
            data = self.parser.parse(self.path)
            self.cache.store(self.id, data)
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
        raise ValueError(f"Unknown parser type: {parser_type}")


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
            library_file = get_data_file("library.yaml")
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
