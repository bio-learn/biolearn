import os
import hashlib
import requests
from urllib.parse import urlparse
import shutil
import appdirs
import pandas as pd
from importlib.resources import open_text


def getOlinkManifest():
    """
    Loads the Olink protein manifest from the embedded data directory.

    Returns:
        pd.DataFrame: A DataFrame containing the Olink protein manifest with UniProt IDs and gene names.
    """
    # Import manifest of Olink proteins
    file = open_text("biolearn.data", "Olink_uniprot_IDs.csv")
    data = pd.read_csv(file)
    # Rename genes
    data["Gene"] = data["Gene"].str.replace("ERVV-1", "ERVV_1")
    data["Gene"] = data["Gene"].str.replace("HLA-", "HLA_")
    return data


# Cache for mappings
_mapping_dict = None
_remap_dict = None


def load_mappings():
    """
    Loads the UniProt ID to gene name mappings from the embedded CSV file.

    Returns:
        tuple: A tuple containing two dictionaries:
            - mapping_dict: Maps UniProt IDs to gene names.
            - remap_dict: Maps gene names to UniProt IDs.
    """
    global _mapping_dict, _remap_dict
    if _mapping_dict is None or _remap_dict is None:
        # Read the file
        file = open_text("biolearn.data", "Olink_uniprot_IDs.csv")
        uniprot_manifest = pd.read_csv(file)

        # Rename genes
        uniprot_manifest["Gene"] = uniprot_manifest["Gene"].str.replace(
            "NT-proBNP", "NTproBNP"
        )
        uniprot_manifest["Gene"] = uniprot_manifest["Gene"].str.replace(
            "ERVV-1", "ERVV_1"
        )
        uniprot_manifest["Gene"] = uniprot_manifest["Gene"].str.replace(
            "HLA-", "HLA_"
        )

        # Create dictionaries
        _mapping_dict = dict(
            zip(uniprot_manifest["UniProt ID"], uniprot_manifest["Gene"])
        )
        _remap_dict = dict(
            zip(uniprot_manifest["Gene"], uniprot_manifest["UniProt ID"])
        )

    return _mapping_dict, _remap_dict


# Map UniProt IDs to gene names
def uniprot_to_gene(data):
    """
    Maps UniProt IDs in the DataFrame to gene names using the preloaded mappings.

    Args:
        data (pd.DataFrame): DataFrame containing UniProt IDs as columns.

    Returns:
        pd.DataFrame: DataFrame with UniProt IDs replaced by gene names.
    """
    mapping_dict, _ = load_mappings()
    data.rename(columns=mapping_dict, inplace=True)
    return data


# Map gene names to UniProt IDs
def gene_to_uniprot(data):
    """
    Maps gene names in the DataFrame to UniProt IDs using the preloaded mappings.

    Args:
        data (pd.DataFrame): DataFrame containing gene names as columns.

    Returns:
        pd.DataFrame: DataFrame with gene names replaced by UniProt IDs.
    """
    _, remap_dict = load_mappings()
    data.rename(columns=remap_dict, inplace=True)
    return data


def get_data_file(relative_path):
    """
    Constructs the full path to a data file located in the biolearn embedded data directory

    Args:
        relative_path (str): The relative path to the data file within the 'data' directory.

    Returns:
        str: The full file path to the specified data file.
    """
    script_dir = os.path.dirname(__file__)
    data_file_path = os.path.join(script_dir, "data", relative_path)
    return data_file_path


def get_test_data_file(relative_path):
    """
    Constructs the full path to a data file located in the biolearn embedded test data directory

    Args:
        relative_path (str): The relative path to the test data file within the 'test/data' directory.

    Returns:
        str: The full file path to the specified test data file.
    """
    script_dir = os.path.dirname(__file__)
    data_file_path = os.path.join(script_dir, "test/data", relative_path)
    return data_file_path


def load_test_data_file(relative_path):
    """
    Loads a test data file from the biolearn embedded test data directory as a dataframe

    Args:
        relative_path (str): The relative path to the test data file within the 'test/data' directory.

    Returns:
        pd.DataFrame: The loaded test data as a DataFrame.
    """
    test_sample = pd.read_csv(get_test_data_file(relative_path), index_col=0)
    return test_sample


def is_url(s):
    """
    Checks if a given string is a valid URL.

    Args:
        s (str): The string to check.

    Returns:
        bool: True if the string is a URL, False otherwise.
    """
    parsed = urlparse(s)
    return bool(parsed.scheme)


def cached_download(url_or_filepath):
    """
    Downloads a file from a URL if not already cached, or verifies the existence of a local file.

    For URLs, the file is downloaded and saved in the user's cache directory using a hash of the URL.
    If the input is a local file path, it checks if the file exists. If not, raises a FileNotFoundError.

    Args:
        url_or_filepath (str): The URL or local file path to download or locate.

    Returns:
        str: The path to the downloaded or located file.

    Raises:
        FileNotFoundError: If the input is a local path and the file does not exist.
        Exception: If the file download fails due to network or HTTP errors.
    """
    if is_url(url_or_filepath):
        url = url_or_filepath
        url_path = urlparse(url).path
        ext = os.path.splitext(url_path)[1]
        filename = hashlib.sha256(url.encode()).hexdigest() + ext

        app_name = "bio-learn"
        download_path = appdirs.user_cache_dir(app_name)
        os.makedirs(download_path, exist_ok=True)
        filepath = os.path.join(download_path, filename)

        if not os.path.exists(filepath):
            # Download the file if it doesn't already exist locally
            try:
                response = requests.get(url, stream=True)
                response.raise_for_status()  # Check for HTTP errors
                with open(filepath, "wb") as out_file:
                    shutil.copyfileobj(response.raw, out_file)
            except requests.RequestException as e:
                raise Exception(
                    f"Failed to download the file from {url}. Error: {e}"
                )
    else:
        if not os.path.exists(url_or_filepath):
            raise FileNotFoundError(
                f"The file does not exist: {url_or_filepath}"
            )
        filepath = url_or_filepath

    return filepath
