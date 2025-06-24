from biolearn.data_library import GeoData


def test_search_returns_df():
    df = GeoData.search(sex="female")
    # DataFrame exists & includes key column
    assert df is not None and not df.empty
    assert "series_id" in df.columns
