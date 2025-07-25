import math
from biolearn.metadata import sex_to_numeric, numeric_to_sex, standardize_sex

def test_sex_numeric_map():
    assert sex_to_numeric("female") == 0
    assert sex_to_numeric("male") == 1
    assert math.isnan(sex_to_numeric("unknown"))
    assert numeric_to_sex(1) == "male"
    assert numeric_to_sex(0) == "female"
    assert numeric_to_sex(float('nan')) == "unknown"

def test_standardize_sex():
    # Test standard cases
    assert standardize_sex("female") == 0
    assert standardize_sex("male") == 1
    assert standardize_sex("f") == 0
    assert standardize_sex("m") == 1
    assert standardize_sex("F") == 0
    assert standardize_sex("M") == 1

    # Test unknown/missing cases
    assert math.isnan(standardize_sex("unknown"))
    assert math.isnan(standardize_sex(""))
    assert math.isnan(standardize_sex(None))
    assert math.isnan(standardize_sex(float('nan')))

    # Test numeric inputs
    assert standardize_sex(0) == 0  # female
    assert standardize_sex(1) == 1  # male
    assert standardize_sex(2) == 0  # GEO encoding: 2=female
    assert math.isnan(standardize_sex(99))  # unknown numeric
