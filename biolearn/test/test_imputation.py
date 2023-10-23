import pytest
import pandas as pd
import numpy as np

from biolearn.imputation import (
    impute_from_standard,
    impute_from_average,
    hybrid_impute,
)

df_test = pd.DataFrame(
    {
        "Sample1": [1, np.nan, 3, 4],
        "Sample2": [np.nan, 1, 2, 3],
        "Sample3": [1, 2, np.nan, 4],
        "Sample4": [2, 1, 3, 4],
        "Sample5": [2, np.nan, 2, 3],
        "Sample6": [4, 2, 2, 4],
    },
    index=["cpg1", "cpg2", "cpg3", "cpg4"],
)

cpg_averages_test = pd.Series(
    {"cpg1": 1.5, "cpg2": 2.5, "cpg3": 3.5, "cpg4": 4.5, "cpg5": 5.5}
)

required_cpgs = ["cpg1", "cpg2", "cpg3", "cpg4", "cpg5"]


def no_missing_values(df):
    return df.isna().sum().sum() == 0


# Tests
def test_impute_from_standard():
    df_filled = impute_from_standard(df_test, cpg_averages_test)
    assert no_missing_values(df_filled)
    assert df_filled.loc["cpg2", "Sample1"] == 2.5


def test_impute_from_standard_specific_cpgs():
    specific_cpgs = ["cpg1", "cpg3"]
    df_filled = impute_from_standard(
        df_test, cpg_averages_test, cpgs_to_impute=specific_cpgs
    )

    print(df_filled)
    print(df_test)
    assert not df_filled.loc[specific_cpgs].isna().any().any()
    assert (
        df_filled.drop(specific_cpgs).isna().sum().sum()
        == df_test.drop(specific_cpgs).isna().sum().sum()
    )
    assert df_filled.loc["cpg1", "Sample2"] == 1.5
    assert df_filled.loc["cpg3", "Sample3"] == 3.5


def test_impute_from_average():
    df_filled = impute_from_average(df_test)
    assert no_missing_values(df_filled)
    assert (
        df_filled.loc["cpg1", "Sample2"] == 2
    )  # Ensure the value was filled with the row average.


def test_impute_from_average_specific_cpgs():
    specific_cpgs = ["cpg1", "cpg3"]
    df_filled = impute_from_average(df_test, cpgs_to_impute=specific_cpgs)

    print(df_filled)
    print(df_test)
    assert not df_filled.loc[specific_cpgs].isna().any().any()
    assert (
        df_filled.drop(specific_cpgs).isna().sum().sum()
        == df_test.drop(specific_cpgs).isna().sum().sum()
    )
    assert df_filled.loc["cpg1", "Sample2"] == 2
    assert df_filled.loc["cpg3", "Sample3"] == 2.4


def test_replacement_with_mean():
    result = hybrid_impute(
        df_test, cpg_averages_test, ["cpg1", "cpg2", "cpg3", "cpg4", "cpg5"]
    )
    assert no_missing_values(result)
    assert result.loc["cpg1", "Sample1"] == 1
    assert result.loc["cpg1", "Sample2"] == 2


def test_replacement_with_source_below_threshold():
    result = hybrid_impute(
        df_test, cpg_averages_test, ["cpg1", "cpg2", "cpg3", "cpg4", "cpg5"]
    )
    assert no_missing_values(result)
    assert result.loc["cpg2", "Sample1"] == 2.5  # Missing value is replaced
    assert result.loc["cpg2", "Sample2"] == 2.5  # Existing value is replaced


def test_addition_of_missing_values():
    result = hybrid_impute(
        df_test, cpg_averages_test, ["cpg1", "cpg2", "cpg3", "cpg4", "cpg5"]
    )
    assert no_missing_values(result)
    assert "cpg5" in result.index
    assert result.loc["cpg5", "Sample1"] == 5.5
    assert result.loc["cpg5", "Sample2"] == 5.5


def test_error_when_required_cpg_not_available():
    with pytest.raises(ValueError):
        hybrid_impute(
            df_test,
            cpg_averages_test.drop("cpg5"),
            ["cpg1", "cpg2", "cpg3", "cpg4", "cpg5"],
        )


def test_modification_of_threshold():
    result = hybrid_impute(
        df_test, cpg_averages_test, ["cpg1", "cpg2", "cpg3", "cpg4", "cpg5"], 1
    )
    # cpg2 for Sample1 should be replaced with source value because threshold is 1
    assert result.loc["cpg2", "Sample1"] == 2.5
    result = hybrid_impute(
        df_test,
        cpg_averages_test,
        ["cpg1", "cpg2", "cpg3", "cpg4", "cpg5"],
        0.5,
    )
    # cpg2 for Sample1 should be replaced with mean value because threshold is 0.5
    assert (
        result.loc["cpg2", "Sample1"] == 1.5
    )  # mean of available value and source value for cpg2


# Run the tests
if __name__ == "__main__":
    pytest.main([__file__])
