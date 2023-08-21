import yaml
import pandas as pd
from biolearn.util import cached_dowload, get_data_file


class GeoMatrixParser:
    def __init__(self, data):
        self.id_row = data.get("id-row")
        self.age_row = data.get("age-row")
        self.age_parse = data.get("age-parse")
        self.matrix_start = data.get("matrix-start")

    def parse(self, file_path):
        ages = pd.read_table(
            file_path,
            index_col=0,
            skiprows=lambda x: x != self.age_row - 1 and x != self.id_row - 1,
        ).transpose()

        dnam = pd.read_table(
            file_path, index_col=0, skiprows=self.matrix_start - 1
        ).transpose()
        dnam["age"] = ages["!Sample_characteristics_ch1"].str[-2:].astype(int)
        dnam = dnam.drop(["!series_matrix_table_end"], axis=1)
        dnam.index.name = "id"
        return dnam


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
