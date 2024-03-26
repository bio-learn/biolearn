import os
import hashlib
import requests
from urllib.parse import urlparse
import shutil
import appdirs
import pandas as pd


def get_data_file(relative_path):
    script_dir = os.path.dirname(
        __file__
    )  # get the directory of the current script
    data_file_path = os.path.join(
        script_dir, "data", relative_path
    )  # build the path to the data file
    return data_file_path


def get_test_data_file(relative_path):
    script_dir = os.path.dirname(
        __file__
    )  # get the directory of the current script
    data_file_path = os.path.join(
        script_dir, "test/data", relative_path
    )  # build the path to the data file
    return data_file_path


def load_test_data_file(relative_path):
    test_sample = pd.read_csv(get_test_data_file(relative_path), index_col=0)
    return test_sample


def is_url(s):
    """Check if the input string is a URL."""
    parsed = urlparse(s)
    return bool(parsed.scheme)


def cached_download(url_or_filepath):
    """Downloads the file at a URL and saves it locally. If called again with the same URL, it will use the saved file.
    If the input is a local file path and the file does not exist, raises a FileNotFoundError.
    Returns the local filepath."""
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
