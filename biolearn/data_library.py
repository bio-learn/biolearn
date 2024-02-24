import yaml
import pandas as pd
import numpy as np
import re
from biolearn.util import cached_download, get_data_file
from IPython.display import display, HTML


def parse_after_colon(s):
    """Extract and return the substring after the first colon."""

    if isinstance(s, str) and s.strip() and ":" in s:
        return s.split(":")[1].strip()

    return ""


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

    def __init__(self, metadata, dnam):
        """
        Initializes the GeoData instance.

        Args:
            metadata (DataFrame): Metadata associated with genomic samples.
            dnam (DataFrame): Methylation data associated with genomic samples.
        """
        self.metadata = metadata
        self.dnam = dnam

    def copy(self):
        """
        Creates a deep copy of the GeoData instance.

        Returns:
            GeoData: A new instance of GeoData with copies of the metadata and dnam DataFrames.
        """
        return GeoData(
            self.metadata.copy(deep=True), self.dnam.copy(deep=True)
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


class GeoMatrixParser:
    parsers = {
        "numeric": lambda s: extract_numeric(parse_after_colon(s)),
        "string": lambda s: parse_after_colon(s),
        "sex": lambda s: sex_parser(parse_after_colon(s)),
    }
    seperators = {"space": " ", "comma": ","}

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

    def parse(self, file_path):
        load_list = self._metadata_load_list()
        load_rows = [x[1] for x in load_list]
        column_names = [x[0] for x in load_list]
        metadata = pd.read_table(
            file_path,
            index_col=0,
            skiprows=lambda x: x != self.id_row - 1 and x not in load_rows,
        )
        metadata.index = column_names
        metadata = metadata.transpose()
        metadata.index.name = "id"
        for col in metadata.columns:
            parser_name = self.metadata[col]["parse"]
            parser = self.parsers[parser_name]
            metadata[col] = metadata[col].apply(parser)
        if self.matrix_start:
            dnam = pd.read_table(
                file_path, index_col=0, skiprows=self.matrix_start - 1
            )
            dnam = dnam.drop(["!series_matrix_table_end"], axis=0)
            dnam.index.name = "id"
        elif self.matrix_file:
            matrix_file_path = cached_download(self.matrix_file)
            print(
                f"Note: This dataset will take a few minutes to load and may fail if you have insufficient memory"
            )
            df = pd.read_csv(
                matrix_file_path, index_col=0, sep=self.matrix_file_seperator
            )
            methylation_df = df.iloc[:, ::2]
            pval_df = df.iloc[:, 1::2]
            pval_df = pval_df.replace("<1E-16", "0", regex=False).astype(
                float, errors="ignore"
            )
            pval_df = pval_df.map(lambda x: np.nan if x > 0.05 else 0)
            # NaN values in pval_df will cause corresponding values in methylation_df to be NaN
            dnam = methylation_df + pval_df.values
            dnam = self._remap_and_prune_columns(dnam, file_path)
        return GeoData(metadata, dnam)

    def _remap_and_prune_columns(self, data, matrix_file_path):
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

    def __init__(self, data):
        """
        Initializes the DataSource instance. Intended to be initialized with data from a library YAML file.

        Args:
            data (dict): A dictionary containing the data source's properties.

        Raises:
            ValueError: If any of the required fields ('id', 'path', 'parser') are missing.
        """
        for field, error_message in self.REQUIRED_FIELDS.items():
            setattr(self, field, data.get(field))
            if getattr(self, field) is None:
                raise ValueError(error_message)

        self.title = data.get("title")
        self.summary = data.get("summary")
        self.format = data.get("format")
        self.organism = data.get("organism")

        self.parser = self._create_parser(data["parser"])

    def load(self):
        """
        Loads the data from the source.

        Returns:
            GeoData: An instance of the GeoData class containing the parsed geographical data.
        """
        file_path = cached_download(self.path)
        return self.parser.parse(file_path)

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
        if parser_type == "geo-matrix":
            return GeoMatrixParser(parser_data)
        # Add more parsers as needed
        raise ValueError(f"Unknown parser type: {parser_type}")


def parse_library_file(library_file):
    data = yaml.safe_load(library_file)
    if data is None:
        raise ValueError("File must be YAML format with 'items' at root")
    if "items" in data:
        data_sources = [DataSource(item) for item in data["items"]]
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

    def __init__(self, library_file=None):
        """
        Initializes the DataLibrary instance.

        Args:
            library_file (str, optional): The file path of the library file. If None,
                                          the biolearn default library will be loaded.

        """
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
            data_sources = parse_library_file(f)
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
