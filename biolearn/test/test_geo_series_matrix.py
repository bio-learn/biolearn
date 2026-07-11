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


from biolearn.data_library import load_geo_metadata


def test_metadata_by_key_ignores_line_numbers():
    series = _series()
    filekey = {"age": {"key": "age", "parse": "numeric"}}
    meta = load_geo_metadata(series, filekey, id_row=33)
    assert list(meta.index) == series.sample_ids()
    assert meta["age"].notna().all()


def test_metadata_by_tag():
    series = _series()
    filekey = {"title": {"tag": "!Sample_title", "parse": "string"}}
    meta = load_geo_metadata(series, filekey, id_row=33)
    assert meta["title"].notna().all()


def test_metadata_offset_corrects_stale_rows():
    series = _series()
    # age characteristic is on line 47. Pretend the config was written when
    # the header sat two lines higher: id_row 31 (real is 33) and age row 45.
    filekey = {"age": {"row": 45, "parse": "numeric"}}
    meta = load_geo_metadata(series, filekey, id_row=31)
    assert meta["age"].notna().all()


def test_metadata_all_unparseable_raises_when_corrected():
    series = _series()
    # With a +2 offset applied, a bad row that lands on a digit-free
    # characteristic must fail loudly rather than return a NaN column.
    # row 40 + offset 2 = line 42, "sample type: whole blood".
    filekey = {"age": {"row": 40, "parse": "numeric"}}
    with pytest.raises(ValueError):
        load_geo_metadata(series, filekey, id_row=31)


def test_legacy_offset_zero_does_not_validate():
    series = _series()
    # Offset 0, legacy row that happens to be non-numeric stays permissive
    # so the 40 unchanged datasets keep their exact behavior.
    filekey = {"plate": {"row": 43, "parse": "string"}}
    meta = load_geo_metadata(series, filekey, id_row=33)
    assert "plate" in meta.columns


from biolearn.data_library import GeoMatrixParser


def test_geo_matrix_parser_id_row_optional():
    # id-row is now optional; auto-detected from the series matrix
    parser = GeoMatrixParser({"type": "geo-matrix", "matrix-start": 74})
    assert parser.id_row is None


from biolearn.data_library import ChallengeDataParser


def test_challenge_parser_accepts_key_tag_and_optional_id_row():
    parser = ChallengeDataParser(
        {
            "matrix-file": "ftp://example/betas.csv.gz",
            "matrix-file-key-tag": "!Sample_description",
            "metadata": {},
        }
    )
    assert parser.matrix_file_key_tag == "!Sample_description"
    assert parser.id_row is None
