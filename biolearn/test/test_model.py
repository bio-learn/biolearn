import pytest
from math import isclose
import pandas as pd
import numpy as np
from biolearn import model
from biolearn.util import (
    get_test_data_file,
    load_test_data_file,
    get_data_file,
)
import pickle


sample_results = load_test_data_file("expected_clock_output.csv")
sample_inputs = load_test_data_file("external/DNAmTestSet.csv")


@pytest.mark.parametrize(
    "model_name, model_entry", model.model_definitions.items()
)
def test_models(model_name, model_entry):
    test_model = model.LinearMethylationModel.from_definition(model_entry)
    actual_results = test_model.predict(sample_inputs).sort_index()

    expected_results = sample_results[model_name].sort_index()

    assert len(expected_results) == len(
        actual_results
    ), f"For {model_name}: DataFrames do not have the same length"

    discrepancies = [
        (idx, expected_results.loc[idx], actual_results.loc[idx])
        for idx in expected_results.index
        if not isclose(
            actual_results.loc[idx], expected_results.loc[idx], abs_tol=1e-5
        )
    ]

    # Check if any discrepancies were found and output them
    assert not discrepancies, "\n".join(
        [
            f"For {model_name}: Discrepancy at index {idx}: expected {expected_age}, but got {actual_age}"
            for idx, expected_age, actual_age in discrepancies
        ]
    )


def test_dunedin_pace_normalization():
    actual = model.dunedin_pace_normalization(sample_inputs)
    data_file_path = get_test_data_file("pace_normalized.pkl")
    with open(data_file_path, "rb") as file:
        expected = pickle.load(file)

    # Finding mismatches based on tolerance
    mask = np.abs(actual - expected) > 0.00000001
    mismatches = actual[mask].stack()

    total_mismatches = mismatches.size
    percentage_mismatched = (total_mismatches / actual.size) * 100

    # Display the mismatches
    for idx, value in enumerate(mismatches.items()):
        if idx == 100:
            break
        print(
            f"Location: {value[0]}, Actual: {actual.at[value[0]]}, Expected: {expected.at[value[0]]}"
        )

    print(
        f"Total mismatches: {total_mismatches} ({percentage_mismatched:.2f}%)"
    )

    # Your actual assertion can be here based on your needs, e.g.,
    assert (
        total_mismatches == 0
    ), "Dataframes are not equal within the given tolerance."


# Run the test
if __name__ == "__main__":
    pytest.main([__file__])
