import pandas as pd
import pytest
from unittest.mock import Mock, patch

from biolearn.data_library import (
    AutoScanGeoMatrixParser,
    DataLibrary,
    NoMatrixDataError,
)
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

    # geo2r now answers scripted requests with a CAPTCHA page, so stub the
    # metadata call and exercise only the matrix check under test
    stub_metadata = pd.DataFrame({"id": ["GSE121633"]})
    with patch.object(
        AutoScanGeoMatrixParser, "_create_metadata", return_value=stub_metadata
    ):
        with pytest.raises(NoMatrixDataError):
            library.get("GSE121633").load()


def test_metadata_query_non_json_response_error():
    parser = AutoScanGeoMatrixParser(
        {
            "matrix_file": "unused",
            "metadata_keys_parse": {},
            "metadata_query": "unused",
        }
    )
    html_response = Mock()
    html_response.json.side_effect = ValueError("Expecting value")

    with patch(
        "biolearn.data_library.requests.get", return_value=html_response
    ):
        with pytest.raises(ValueError, match="did not return JSON"):
            parser._create_metadata("https://example.org/geo2r")
