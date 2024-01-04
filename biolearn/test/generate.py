from biolearn.data_library import DataSource
import pandas as pd
import os
import sys
import hashlib

source_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41169/matrix/GSE41169_series_matrix.txt.gz"

# Paths for the script, data, and metadata files
script_dir = os.path.dirname(__file__)
data_folder_path = os.path.join(script_dir, "data/external/")
data_file_path = os.path.join(data_folder_path, "DNAmTestSet.csv")
metadata_file_path = os.path.join(data_folder_path, "testset_metadata.csv")


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
    print("Loading data source")
    data_source_spec = {
        "id": "TestData",
        "path": source_url,
        "parser": {
            "type": "geo-matrix",
            "id-row": 33,
            "metadata": {
                "age": {"row": 47, "parse": "numeric"},
                "sex": {"row": 41, "parse": "sex"},
            },
            "matrix-start": 73,
        },
    }
    source = DataSource(data_source_spec)
    loaded_data = source.load()

    print("Processing DNAm data")
    dnam = (
        loaded_data.dnam.transpose().head(10).transpose()
    )  # Transpose and limit to 10 items
    dnam.to_csv(data_file_path)

    print("Verifying new DNAm data file")
    if not is_file_valid(data_file_path, expected_hash):
        print(
            "ERROR: Generated DNAm data file does not match expected. Exiting."
        )
        sys.exit(1)
    else:
        print("DNAm data file matches existing hash. Done.")

    print("Processing metadata")
    metadata = loaded_data.metadata.head(10)  # Limit to 10 items
    metadata.to_csv(metadata_file_path)
    print("Metadata file generated and saved.")
else:
    print("File matches existing hash. Done.")
