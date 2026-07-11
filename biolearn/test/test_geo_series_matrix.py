from biolearn.data_library import GeoSeriesMatrix
from biolearn.util import get_test_data_file


def _series():
    return GeoSeriesMatrix(get_test_data_file("geo_dnam_test_file"))


def test_sample_ids_from_geo_accession():
    series = _series()
    ids = series.sample_ids()
    assert ids == [
        "GSM1009660",
        "GSM1009661",
        "GSM1009662",
        "GSM1009663",
        "GSM1009664",
    ]


def test_id_row_is_geo_accession_line():
    assert _series().id_row == 33


def test_matrix_start_is_after_table_begin():
    # table begin marker is on line 73, header row follows on 74
    assert _series().matrix_start == 74


def test_tag_values_returns_sample_tag():
    titles = _series().tag_values("!Sample_title")
    assert len(titles) == 5
    assert _series().tag_values("!Sample_not_present") is None


def test_characteristic_by_key_is_case_insensitive():
    series = _series()
    gender = series.characteristic_values("Gender")
    assert len(gender) == 5
    assert all(v.lower().startswith("gender:") for v in gender)
    assert series.characteristic_values("no-such-key") is None


def test_line_values_is_one_based():
    # line 33 is the geo accession row
    assert _series().line_values(33)[0] == "GSM1009660"


def test_id_offset():
    series = _series()
    assert series.id_offset(33) == 0
    assert series.id_offset(31) == 2
    assert series.id_offset(None) == 0


import pytest
import pandas as pd
from biolearn.data_library import build_column_mapping, map_and_prune_columns


def test_build_column_mapping_by_tag():
    series = _series()
    mapping = build_column_mapping(series, key_tag="!Sample_title")
    assert len(mapping) == 5
    assert all(str(v).startswith("GSM") for v in mapping.values())


def test_build_column_mapping_by_offset_corrected_line():
    series = _series()
    # !Sample_title is line 32; feed a stale line 30 with a +2 offset
    mapping = build_column_mapping(series, key_line=30, offset=2)
    assert len(mapping) == 5
    assert all(str(v).startswith("GSM") for v in mapping.values())


def test_build_column_mapping_missing_source_raises():
    with pytest.raises(ValueError):
        build_column_mapping(_series(), key_tag="!Sample_not_present")


def test_map_and_prune_keeps_only_mapped_columns():
    data = pd.DataFrame(
        {"sentrixA": [1, 2], "sentrixB": [3, 4], "junk": [5, 6]}
    )
    mapping = {"sentrixA": "GSM1", "sentrixB": "GSM2"}
    pruned = map_and_prune_columns(data, mapping)
    assert list(pruned.columns) == ["GSM1", "GSM2"]
