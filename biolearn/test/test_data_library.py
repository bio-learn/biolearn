import pytest
import os
from io import StringIO
from biolearn.load import parse_library_file

def load_test_data_file(relative_path):
    script_dir = os.path.dirname(
        __file__
    )  # get the directory of the current script
    data_file_path = os.path.join(
        script_dir, "data", relative_path
    )  # build the path to the data file
    test_sample = pd.read_csv(data_file_path, index_col=0)
    return test_sample

def test_blank_file_gives_error():
    file_contents = " "
    test_file = StringIO()
    test_file.write(file_contents)
    with pytest.raises(ValueError) as e:
        parse_library_file(test_file)
    assert str(e.value) == "File must be YAML format with items at root"

