import pytest
from math import isclose
import pandas as pd
from collections import defaultdict
import numpy as np
from biolearn import model
from biolearn.util import (
    get_test_data_file,
    load_test_data_file,
)
from biolearn.data_library import GeoData
import pickle

TOLERANCES = defaultdict(lambda: 1e-5)
# AltumAge: TF→Torch port can differ ~1e-4 across platforms/versions.
TOLERANCES["AltumAge"] = 2e-4
# MiAge: iterative algorithm can drift ~3e-5 across Python/platform versions.
TOLERANCES["MiAge"] = 5e-5

sample_inputs = load_test_data_file("testset/testset_methylation_part0.csv")


@pytest.mark.parametrize(
    "model_name, model_entry", model.model_definitions.items()
)
def test_models(model_name, model_entry):
    # TODO: Make testing more intelligent for different model types.
    # TODO: Add testing for LinearTranscriptomicModel
    # Skip models that don't have tests
    model_type = model_entry["model"]["type"]
    if model_type in [
        "NotImplemented",
        "LinearTranscriptomicModel",
        "HurdleAPIModel",
    ]:
        pytest.skip(
            f"Model type {model_type} for {model_name} does not have a testing pattern - skipping test"
        )

    test_data = GeoData.load_csv(
        get_test_data_file("testset/"), "testset", validate=False
    )

    # Check if the model class exists
    try:
        model_class = getattr(model, model_entry["model"]["type"])
    except AttributeError:
        pytest.fail(
            f"Model class '{model_entry['model']['type']}' not found in biolearn.model"
        )

    # Instantiate the model
    test_model = model_class.from_definition(model_entry)

    if model_type == "LinearMethylationModel":
        required_cpgs = test_model.methylation_sites()
        missing_cpgs = sorted(set(required_cpgs) - set(test_data.dnam.index))
        if missing_cpgs:
            missing_preview = missing_cpgs[0]
            with pytest.raises(
                ValueError,
                match=rf"Missing required CpG sites.*{missing_preview}",
            ):
                test_model.predict(test_data)
            return

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
                actual_val, expected_val, abs_tol=TOLERANCES[model_name]
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

    # The reference fixture stores probe IDs in an ``ID_REF`` column with a
    # default RangeIndex, whereas the loader returns CpG-indexed data. Align the
    # fixture's orientation to ``actual`` before comparing; otherwise the axes
    # don't match and the comparison is vacuous. (On older pandas ``stack()``
    # silently dropped the misaligned NaNs, so this test passed without ever
    # comparing values; pandas 3.0 keeps them, surfacing a KeyError instead.)
    if "ID_REF" in expected.columns:
        expected = expected.set_index("ID_REF")
        expected.index.name = actual.index.name
    expected = expected.reindex(index=actual.index, columns=actual.columns)

    # Count mismatches from the boolean mask directly, independent of how
    # ``stack()`` treats NaNs across pandas versions.
    mismatch_mask = np.abs(actual - expected) > 0.000001
    total_mismatches = int(mismatch_mask.to_numpy().sum())
    percentage_mismatched = (total_mismatches / actual.size) * 100

    # Display up to the first 100 mismatching cells, if any.
    if total_mismatches:
        stacked = mismatch_mask.stack()
        for row, col in stacked[stacked].index[:100]:
            print(
                f"Location: ({row}, {col}), "
                f"Actual: {actual.at[row, col]}, Expected: {expected.at[row, col]}"
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
