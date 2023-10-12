import yaml
import pandas as pd
import re
from biolearn.util import cached_dowload, get_data_file


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
    def __init__(self, metadata, dnam):
        # Metadata should have rows being samples and columns being data fields
        self.metadata = metadata
        # Methylation data should have columns as samples and rows as methylation sites
        self.dnam = dnam


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
    REQUIRED_FIELDS = {
        "id": "'id' key is missing in item",
        "path": "'path' key is missing in item",
        "parser": "'parser' key is missing in item",
    }

    def __init__(self, data):
        for field, error_message in self.REQUIRED_FIELDS.items():
            setattr(self, field, data.get(field))
            if getattr(self, field) is None:
                raise ValueError(error_message)

        self.title = data.get("title")
        self.summary = data.get("summary")
        self.format = data.get("format")
        self.organism = data.get("organism")

        self.parser = self._create_parser(data["parser"])

    def _create_parser(self, parser_data):
        parser_type = parser_data.get("type")
        if parser_type is None:
            raise ValueError("Parser type is missing")
        if parser_type == "geo-matrix":
            return GeoMatrixParser(parser_data)
        # Add more parsers as needed
        raise ValueError(f"Unknown parser type: {parser_type}")

    def load(self):
        file_path = cached_dowload(self.path)
        return self.parser.parse(file_path)

    def __repr__(self):
        return f"DataSource(ID: {self.id}\nTitle: {self.title})"


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
    def __init__(self, library_file=None):
        self.sources = []
        if library_file is None:
            library_file = get_data_file("library.yaml")
        self.load_sources(library_file)

    def load_sources(self, library_file):
        with open(library_file, "r") as f:
            data_sources = parse_library_file(f)
            self.sources.extend(data_sources)

    def get(self, source_id):
        for source in self.sources:
            if source.id == source_id:
                return source
        return None

    def lookup_sources(self, organism=None, format=None):
        matches = []
        for source in self.sources:
            if (organism is None or source.organism == organism) and (
                format is None or source.format == format
            ):
                matches.append(source)
        return matches
