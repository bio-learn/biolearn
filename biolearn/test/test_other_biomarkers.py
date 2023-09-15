from biolearn.util import load_test_data_file
from biolearn.other_biomarkers import estimate_sex
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest


sample_inputs = load_test_data_file("external/DNAmTestSet.csv")


def test_estimate_sex():
    expected = load_test_data_file("expected_estimate_sex_output.csv")
    actual = estimate_sex(sample_inputs)
    print(actual)
    assert_frame_equal(expected, actual)
