from biolearn.data_library import DataSource
import pandas as pd
import os
import sys
import hashlib

source_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41169/matrix/GSE41169_series_matrix.txt.gz"

script_dir = os.path.dirname(
    __file__
)  # get the directory of the current script
data_folder_path = os.path.join(
    script_dir, "data/external/"
)  # build the path to the data file
data_file_path = os.path.join(
    data_folder_path, "DNAmTestSet.csv"
)  # build the path to the data file

expected_hash = (
    "3ee979ae6d2f1572f72f02e951073af0910b5b1576d46177cbfb76932fda6693"
)


def compute_file_hash(file_path, chunk_size=1024 * 1024):
    hash_algorithm = hashlib.sha256()
    with open(file_path, "rb") as f:
        while True:
            data = f.read(chunk_size)
            if not data:
                break
            hash_algorithm.update(data)
    return hash_algorithm.hexdigest()


def is_file_valid(file_path, expected_hash):
    if os.path.exists(file_path):
        actual_hash = compute_file_hash(file_path)
        if actual_hash != expected_hash:
            print(
                f"Warning: Actual file hash {actual_hash} does not equal expected hash {expected_hash}"
            )
        return actual_hash == expected_hash
    return False


def ensure_folder_exists(path):
    # Get the directory part of the path (if it is a file)
    dir_path = os.path.dirname(path)

    if os.path.isdir(path):
        dir_path = path

    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


ensure_folder_exists(data_file_path)

print("Checking if file generation is needed")
if not is_file_valid(data_file_path, expected_hash):
    print("Generating file")
    data_source_spec = {
        "id": "TestData",
        "path": source_url,
        "parser": {
            "type": "geo-matrix",
            "id-row": 33,
            "metadata": {
                "age": {"row": 47, "parse": "numeric"},
            },
            "matrix-start": 73,
        },
    }
    source = DataSource(data_source_spec)
    source_data = source.load().dnam
    source_data.transpose().head(10).transpose().to_csv(data_file_path)
    print("Verifying new file")
    if not is_file_valid(data_file_path, expected_hash):
        print("ERROR: Generated file does not match expected. Exiting.")
        sys.exit(1)
    else:
        print("File matches existing hash. Done.")
else:
    print("File matches existing hash. Done.")
