import pytest

from biolearn.data_library import DataLibrary, NoMatrixDataError
from biolearn.util import get_data_file


def test_can_load_autoscan_library_file():
    library = DataLibrary(
        library_file=get_data_file("geo_autoscan_library.yaml")
    )

    assert len(library.sources) > 1


# TODO: Figure out why this is failing on Github with 403
# def test_series_has_matrix_data():
#     library = DataLibrary(
#         library_file=get_data_file("geo_autoscan_library.yaml")
#     )

#     data = library.get("GSE100386").load()

#     assert len(data.metadata) == len(data.dnam.columns)


def test_series_has_no_matrix_data_error():
    library = DataLibrary(
        library_file=get_data_file("geo_autoscan_library.yaml")
    )

    with pytest.raises(NoMatrixDataError):
        library.get("GSE121633").load()
