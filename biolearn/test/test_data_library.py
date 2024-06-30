import pytest
import os
import pandas as pd
import numpy as np
from io import StringIO
from biolearn.data_library import (
    parse_library_file,
    DataSource,
    DataLibrary,
    GeoData,
)
from biolearn.model import get_data_file
from biolearn.util import get_test_data_file
from biolearn.util import load_test_data_file, get_test_data_file

import pandas as pd


def test_quality_report():
    sample_inputs = load_test_data_file("external/DNAmTestSet.csv")
    sample_metadata = load_test_data_file("external/testset_metadata.csv")
    geo_data_instance = GeoData(sample_metadata, sample_inputs)
    actual_report = geo_data_instance.quality_report()

    expected_report = load_test_data_file("test_quality_report.csv")

    # Directly compare the actual and expected reports
    pd.testing.assert_frame_equal(
        actual_report.sample_report,
        expected_report,
        check_dtype=True,
        check_like=True,
    )


def test_blank_file_gives_error():
    file_contents = " "
    test_file = StringIO()
    test_file.write(file_contents)
    with pytest.raises(ValueError) as e:
        parse_library_file(test_file)
    assert str(e.value) == "File must be YAML format with 'items' at root"


def test_can_load_single_item_correctly():
    file_contents = """
---
items:
- id: GSE40279
  title: An informative title
  summary: A Useful Summary
  format: Illumina450k
  organism: human
  path: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz
  parser:
    type: geo-matrix
    id-row: 12
    metadata:
      age:
        row: 44
        parse: numeric
    matrix-start: 72
"""
    test_file = StringIO(file_contents)
    actual = parse_library_file(test_file)

    assert len(actual) == 1
    item = actual[0]
    assert item.id == "GSE40279"
    assert item.title == "An informative title"
    assert item.format == "Illumina450k"
    assert item.organism == "human"
    assert (
        item.path
        == "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz"
    )
    parser = item.parser
    assert parser.id_row == 12
    assert parser.metadata == {"age": {"row": 44, "parse": "numeric"}}
    assert parser.matrix_start == 72


def test_missing_id_gives_error():
    file_contents = """
---
items:
- title: An informative title
  path: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz
  parser:
    type: geo-matrix
    id-row: 12
    metadata:
      age:
        row: 44
        parse: numeric
    matrix-start: 72
"""
    test_file = StringIO(file_contents)
    with pytest.raises(ValueError) as e:
        parse_library_file(test_file)
    assert str(e.value) == "'id' key is missing in item"


def test_unknown_parser_type_gives_error():
    file_contents = """
---
items:
- id: GSE40279
  path: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz
  parser:
    type: unknown-parser
    id-row: 12
    metadata:
      age:
        row: 44
        parse: numeric
    matrix-start: 72
"""
    test_file = StringIO(file_contents)
    with pytest.raises(ValueError) as e:
        parse_library_file(test_file)
    assert str(e.value) == "Unknown parser type: unknown-parser"


def test_missing_path_gives_error():
    file_contents = """
---
items:
- id: GSE40279
  parser:
    type: geo-matrix
    id-row: 12
    metadata:
      age:
        row: 44
        parse: numeric
    matrix-start: 72
"""
    test_file = StringIO(file_contents)
    with pytest.raises(ValueError) as e:
        parse_library_file(test_file)
    assert str(e.value) == "'path' key is missing in item"


def test_parser_missing_id_row_gives_error():
    file_contents = """
---
items:
- id: GSE40279
  path: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz
  parser:
    type: geo-matrix
    misspelled-id-row: 12
    metadata:
      age:
        row: 44
        parse: numeric
    matrix-start: 72
"""
    test_file = StringIO(file_contents)
    with pytest.raises(ValueError) as e:
        parse_library_file(test_file)
    assert str(e.value) == "Parser not valid: missing id-row"


def test_missing_parser_gives_error():
    file_contents = """
---
items:
- id: GSE40279
  path: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz
"""
    test_file = StringIO(file_contents)
    with pytest.raises(ValueError) as e:
        parse_library_file(test_file)
    assert str(e.value) == "'parser' key is missing in item"


def test_can_load_library_file():
    library_file = get_data_file("library.yaml")
    with open(library_file, "r") as the_file:
        actual = parse_library_file(the_file)

    assert len(actual) > 1


def test_can_load_dnam():
    data_source_spec = {
        "id": "TestData",
        "path": get_test_data_file("geo_dnam_test_file"),
        "parser": {
            "type": "geo-matrix",
            "id-row": 33,
            "metadata": {
                "age": {"row": 47, "parse": "numeric"},
                "sex": {"row": 41, "parse": "sex"},
                "cancer": {"row": 50, "parse": "string"},
            },
            "matrix-start": 74,
        },
    }
    source = DataSource(data_source_spec)

    df = source.load()
    # Verify data set is of known size
    dnam = df.dnam
    assert dnam.shape == (37, 5)
    assert df.metadata.shape == (5, 3)
    assert "cancer" in df.metadata.columns.to_list()
    assert np.issubdtype(df.metadata["age"], np.number)
    assert np.issubdtype(df.metadata["sex"], np.number)
    assert (df.metadata["sex"] != 0).all()
    assert all(
        np.issubdtype(dnam[col].dtype, np.number) for col in dnam.columns
    )


def test_load_dnam_with_matrix_file():
    data_source_spec = {
        "id": "TestMatrixFile",
        "path": get_test_data_file("geo_dnam_test_file"),
        "parser": {
            "type": "geo-matrix",
            "id-row": 33,
            "metadata": {
                "age": {"row": 47, "parse": "numeric"},
                "sex": {"row": 41, "parse": "sex"},
                "cancer": {"row": 50, "parse": "string"},
            },
            "matrix-file": get_test_data_file("test_supplementary_matrix.txt"),
            "matrix-file-seperator": "space",
            "matrix-file-key-line": 72,
            "matrix-file-format": "pvalue",
        },
    }
    source = DataSource(data_source_spec)

    df = source.load()
    dnam = df.dnam
    metadata = df.metadata

    assert dnam.shape == (37, 5)
    assert metadata.shape == (5, 3)
    assert "cancer" in metadata.columns.to_list()
    assert np.issubdtype(metadata["age"], np.number)
    assert np.issubdtype(metadata["sex"], np.number)
    assert (metadata["sex"] != 0).all()
    assert all(
        np.issubdtype(dnam[col].dtype, np.number) for col in dnam.columns
    )
    assert dnam.isna().sum().sum() == 8
    assert all(col in metadata.index for col in dnam.columns)


def test_load_sources_append():
    # Initialize DataLibrary
    library = DataLibrary(
        library_file=get_test_data_file("library_files/library.yaml")
    )

    # Check if sources were loaded
    assert len(library.sources) == 2


def test_get_source_by_id():
    # Initialize DataLibrary
    library = DataLibrary(
        library_file=get_test_data_file("library_files/library.yaml")
    )

    # Get source by ID
    source = library.get("GSE40279")
    assert source.id == "GSE40279"
    assert (
        source.title
        == "Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates"
    )


def test_lookup_sources_by_format():
    # Initialize DataLibrary
    library = DataLibrary(
        library_file=get_test_data_file("library_files/library.yaml")
    )

    matches = library.lookup_sources(organism="human", format="Illumina27k")
    assert len(matches) == 1
    assert matches[0].id == "GSE19711"


def test_lookup_sources_by_organism():
    # Initialize DataLibrary
    library = DataLibrary(
        library_file=get_test_data_file("library_files/library.yaml")
    )

    matches = library.lookup_sources(organism="human")
    assert len(matches) == 2
