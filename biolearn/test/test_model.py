import pytest
from math import isclose
import pandas as pd
import numpy as np
from biolearn import model
from biolearn.util import get_test_data_file, load_test_data_file
from biolearn.data_library import GeoData
import pickle


sample_inputs = load_test_data_file("testset/testset_methylation_part0.csv")


@pytest.mark.parametrize(
    "model_name, model_entry", model.model_definitions.items()
)
def test_models(model_name, model_entry):
    # TODO: Make testing more intelligent for different model types.
    # TODO: Add testing for LinearTranscriptomicModel
    # Skip models that don't have tests
    model_type = model_entry["model"]["type"]
    if model_type in ["NotImplemented"]:
        pytest.skip(
            f"Model type {model_type} for {model_name} does not have a testing pattern - skipping test"
        )

    test_data = GeoData.load_csv(get_test_data_file("testset/"), "testset")

    # Check if the model class exists
    try:
        model_class = getattr(model, model_entry["model"]["type"])
    except AttributeError:
        pytest.fail(
            f"Model class '{model_entry['model']['type']}' not found in biolearn.model"
        )

    # Instantiate the model
    test_model = model_class.from_definition(model_entry)

    actual_results = test_model.predict(test_data).sort_index()

    # Load the expected results
    expected_results = load_test_data_file(
        f"expected_model_outputs/{model_name}.csv"
    ).sort_index()

    # Assertions and discrepancy checks
    assert len(expected_results) == len(
        actual_results
    ), f"For {model_name}: DataFrames do not have the same length"

    discrepancies = []
    for idx in expected_results.index:
        for col in expected_results.columns:
            expected_val = expected_results.loc[idx, col]
            actual_val = actual_results.loc[idx, col]

            if isinstance(expected_val, float) and not isclose(
                actual_val, expected_val, abs_tol=1e-5
            ):
                discrepancies.append((idx, col, expected_val, actual_val))
            elif (
                not isinstance(expected_val, float)
                and expected_val != actual_val
            ):
                print(type(expected_val))
                discrepancies.append((idx, col, expected_val, actual_val))

    assert not discrepancies, "\n".join(
        [
            f"For {model_name}: Discrepancy at index {idx}, column {col}: expected {expected}, but got {actual}"
            for idx, col, expected, actual in discrepancies
        ]
    )


def test_dunedin_pace_normalization():

    test_data = GeoData.load_csv(get_test_data_file("testset/"), "testset")
    sample_inputs = test_data.dnam

    actual = model.dunedin_pace_normalization(sample_inputs)
    data_file_path = get_test_data_file("pace_normalized.pkl")
    with open(data_file_path, "rb") as file:
        expected = pickle.load(file)

    # Finding mismatches based on tolerance
    mask = np.abs(actual - expected) > 0.000001
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
