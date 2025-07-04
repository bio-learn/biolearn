from biolearn.data_library import DataLibrary


def test_search_returns_df():
    library = DataLibrary()
    df = library.search(sex="female")
    # DataFrame exists & includes key column
    assert df is not None and not df.empty
    assert "series_id" in df.columns
