from biolearn.load import load_dnam
import pandas as pd
import os
import sys
import hashlib

source_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41169/matrix/GSE41169_series_matrix.txt.gz"

script_dir = os.path.dirname(
    __file__
)  # get the directory of the current script
data_file_path = os.path.join(
    script_dir, "data/external/DNAmTestSet.csv"
)  # build the path to the data file

expected_hash = (
    "62b8b9744c42f6c5bc1f0a40eac4bccfeebc043c7c33da1c0a4eed395bf19003"
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
        return compute_file_hash(file_path) == expected_hash
    return False


print("Checking if file generation is needed")
if not is_file_valid(data_file_path, expected_hash):
    print("Generating file")
    source_data = load_dnam(
        dnam_file=source_url, id_row=32, age_row=46, skiprows=72
    )
    source_data.head(10).transpose().to_csv(data_file_path)
    print("Verifying new file")
    if not is_file_valid(data_file_path, expected_hash):
        print("ERROR: Generated file does not match expected. Exiting.")
        sys.exit(1)
    else:
        print("File matches existing hash. Done.")
else:
    print("File matches existing hash. Done.")
