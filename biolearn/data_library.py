import yaml
import pandas as pd
import re
from biolearn.util import cached_download, get_data_file


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


class GeoMatrixParser:
    parsers = {
        "numeric": lambda s: extract_numeric(parse_after_colon(s)),
        "string": lambda s: parse_after_colon(s),
        "sex": lambda s: sex_parser(parse_after_colon(s)),
    }

    def __init__(self, data):
        if data.get("id-row") is None:
            raise ValueError("Parser not valid: missing id-row")
        self.id_row = data.get("id-row")
        self.metadata = data.get("metadata")
        self.matrix_start = data.get("matrix-start")

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
        dnam = pd.read_table(
            file_path, index_col=0, skiprows=self.matrix_start - 1
        )
        dnam = dnam.drop(["!series_matrix_table_end"], axis=0)
        dnam.index.name = "id"
        return GeoData(metadata, dnam)

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
