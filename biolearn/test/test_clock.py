import pytest
import os
from math import isclose
import pandas as pd
from biolearn import clock


def load_test_data_file(relative_path):
    script_dir = os.path.dirname(
        __file__
    )  # get the directory of the current script
    data_file_path = os.path.join(
        script_dir, "data", relative_path
    )  # build the path to the data file
    test_sample = pd.read_csv(data_file_path, index_col=0)
    return test_sample

sample_results = load_test_data_file("expected_clock_output.csv")
sample_inputs = load_test_data_file("external/DNAmTestSet.csv").transpose()

def check_clock_against_sample(clock_function, results_column_name):
    expected_results = sample_results[results_column_name].sort_index()
    actual_results = clock_function(sample_inputs)['biological_age'].sort_index()
    assert len(expected_results) == len(actual_results), "DataFrames do not have the same length"
    test_passed = True

    for idx in expected_results.index:
        expected_age = expected_results.loc[idx]
        actual_age = actual_results.loc[idx]

        try:
            assert isclose(actual_age, expected_age, abs_tol=1e-2)
        except AssertionError:
            test_passed = False
            print(f"Discrepancy at index {idx}: expected {expected_age}, but got {actual_age}")

    return test_passed

def test_horvathv1_sample():
    assert check_clock_against_sample(clock.horvathv1, "Horvath1")

def test_horvathv2_sample():
    assert check_clock_against_sample(clock.horvathv2, "Horvath2")

def test_hannum_sample():
    assert check_clock_against_sample(clock.hannum, "Hannum")

def test_phenoage_sample():
    assert check_clock_against_sample(clock.phenoage, "PhenoAge")

# Results missing from expected file
# def test_bohlin_sample():
#     assert check_clock_against_sample(clock.bohlin, "Bohlin")

def test_alcohol_mccartney_sample():
    assert check_clock_against_sample(clock.alcohol_mccartney, "Alcohol_McCartney")

def test_bmi_mccartney_sample():
    assert check_clock_against_sample(clock.bmi_mccartney, "BMI_McCartney")

def test_dnam_tl_sample():
    assert check_clock_against_sample(clock.dnam_tl, "DNAmTL")

# Results missing from expected file
# def test_dnam_clock_cortical_sample():
#     assert check_clock_against_sample(clock.dnam_clock_cortical, "DNAmClockCortical")

def test_hrs_in_ch_phenoage():
    assert check_clock_against_sample(clock.hrs_in_ch_phenoage, "HRSInChPhenoAge")

def test_knight():
    assert check_clock_against_sample(clock.knight, "Knight")

def test_lee_control_sample():
    assert check_clock_against_sample(clock.lee_control, "LeeControl")

def test_lee_refined_robust_sample():
    assert check_clock_against_sample(clock.lee_refined_robust, "LeeRefinedRobust")

def test_lee_robust_sample():
    assert check_clock_against_sample(clock.lee_robust, "LeeRobust")

def test_lin_sample():
    assert check_clock_against_sample(clock.lin, "Lin")

# Coeffecient file is broken
# def test_mi_age_sample():
#     assert check_clock_against_sample(clock.mi_age, "MiAge")

def test_pedbe_sample():
    assert check_clock_against_sample(clock.pedbe, "PEDBE")

def test_smoking_mccartney_sample():
    assert check_clock_against_sample(clock.smoking_mccartney, "Smoking_McCartney")

def test_vidal_bralo_sample():
    assert check_clock_against_sample(clock.vidal_bralo, "VidalBralo")

def test_zhang_10_sample():
    assert check_clock_against_sample(clock.zhang_10, "Zhang")

def test_zhang_2019_sample():
    assert check_clock_against_sample(clock.zhang_2019, "Zhang2019")

def test_mayne_sample():
    assert check_clock_against_sample(clock.lee_robust, "Mayne")

# Run the test
if __name__ == "__main__":
    pytest.main([__file__])
