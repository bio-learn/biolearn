import os
import hashlib
import requests
from urllib.parse import urlparse
import shutil
import appdirs
import pandas as pd


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
