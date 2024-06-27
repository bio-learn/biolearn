import pytest

from biolearn.data_library import (
    DataLibrary,
    NoMatrixDataError
)
from biolearn.util import get_test_data_file


def test_can_load_new_yaml_file():
    library = DataLibrary(
        library_file=get_test_data_file("library_files/geo_entry_library.yaml")
    )

    assert len(library.sources) == 574


def test_series_has_matrix_data():
    library = DataLibrary(
        library_file=get_test_data_file("library_files/geo_entry_library.yaml")
    )

    data = library.get("GSE100386").load()

    assert len(data.metadata) == len(data.dnam.columns)


def test_series_has_no_matrix_data_error():
    library = DataLibrary(
        library_file=get_test_data_file("library_files/geo_entry_library.yaml")
    )

    with pytest.raises(NoMatrixDataError):
        library.get("GSE121633").load()
