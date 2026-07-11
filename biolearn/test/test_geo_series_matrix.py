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
