import pytest
import os
import pandas as pd
import numpy as np
from io import StringIO
from biolearn.data_library import parse_library_file, DataSource
from biolearn.clock import get_data_file


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
    age-row: 44
    age-parse: after-colon
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
    assert parser.age_row == 44
    assert parser.age_parse == "after-colon"
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
    age-row: 44
    age-parse: after-colon
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
    age-row: 44
    age-parse: after-colon
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
    age-row: 44
    age-parse: after-colon
    matrix-start: 72
"""
    test_file = StringIO(file_contents)
    with pytest.raises(ValueError) as e:
        parse_library_file(test_file)
    assert str(e.value) == "'path' key is missing in item"


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
    script_dir = os.path.dirname(
        __file__
    )  # get the directory of the current script
    data_file_path = os.path.join(
        script_dir, "data", "geo_dnam_test_file"
    )  # build the path to the data file

    data_source_spec = {
        "id": "TestData",
        "path": data_file_path,
        "parser": {
            "type": "geo-matrix",
            "id-row": 33,
            "age-row": 47,
            "age-parse": "after-colon",
            "matrix-start": 73,
        },
    }
    source = DataSource(data_source_spec)
    
    df = source.load()
    # Verify data set is of known size
    assert df.shape == (5, 38)
    assert "age" in df.columns.to_list()
    assert all(np.issubdtype(df[col].dtype, np.number) for col in df.columns)
