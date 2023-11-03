from biolearn import load
import pytest


load_columns = ["sex", "age", "glucose", "is_dead", "months_until_death"]


def test_fhs_columns():
    df = load.load_fhs()
    verify_expected_columns(df)


def test_nhanes_columns():
    df = load.load_nhanes(2010)
    verify_expected_columns(df)


def test_can_load_nhanes_2012():
    df = load.load_nhanes(2012)


def test_expected_error_when_loading_unsupported_year_nhanes():
    with pytest.raises(ValueError):
        df = load.load_nhanes(1913)


def verify_expected_columns(df):
    actual_columns = set(df.columns.to_list())
    missing_columns = set(load_columns) - actual_columns
    extra_columns = actual_columns - set(load_columns)
    assert (
        len(missing_columns) == 0
    ), f"Missing expected columns: {missing_columns} \n Found extra columns: {extra_columns}"


# Run the test
if __name__ == "__main__":
    pytest.main([__file__])
