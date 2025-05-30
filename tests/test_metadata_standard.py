from biolearn.metadata import sex_to_numeric, numeric_to_sex

def test_sex_numeric_map():
    assert sex_to_numeric("female") == 0
    assert sex_to_numeric("male") == 1
    assert sex_to_numeric("unknown") == -1
    assert numeric_to_sex(1) == "male"
