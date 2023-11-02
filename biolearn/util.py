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


def cached_download(url_or_filepath):
    """Downloads the file at a URL and saves it locally. If called again with the same URL it will use the saved file. Returns the local filepath"""
    # Hash the URL to create a unique filename
    if os.path.isfile(url_or_filepath):
        # If the provided URL is a local file path, return it directly
        return url_or_filepath

    url = url_or_filepath
    url_path = urlparse(url).path
    ext = os.path.splitext(url_path)[1]
    filename = hashlib.sha256(url.encode()).hexdigest() + ext

    app_name = "bio-learn"
    download_path = appdirs.user_cache_dir(app_name)

    # Ensure download path exists
    os.makedirs(download_path, exist_ok=True)

    filepath = os.path.join(download_path, filename)

    if os.path.exists(filepath):
        # If the file is already downloaded, return the file path
        return filepath
    else:
        # Try to download the file
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an HTTPError if one occurred

        # If the file is not downloaded yet, download and save it
        with open(filepath, "wb") as out_file:
            shutil.copyfileobj(response.raw, out_file)

        # Return the file path
        return filepath
