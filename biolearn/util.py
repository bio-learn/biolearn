import hashlib
import os
from importlib.resources import open_text
from urllib.parse import urlparse

import appdirs
import pandas as pd
import requests


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
    if is_url(relative_path):
        # If the input is a URL, download and cache the file
        return cached_download(relative_path)
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


def cached_download(
    url_or_filepath,
    show_progress=False,
    progress_label=None,
    force_download=False,
):
    """
    Downloads a file from a URL if not already cached, or verifies the existence of a local file.

    For URLs, the file is downloaded and saved in the user's cache directory using a hash of the URL.
    If a partial download exists, it attempts to resume with an HTTP Range request.
    If the input is a local file path, it checks if the file exists. If not, raises a FileNotFoundError.

    Args:
        url_or_filepath (str): The URL or local file path to download or locate.
        show_progress (bool): If True, prints periodic download progress.
        progress_label (str | None): Optional label for progress output.
        force_download (bool): If True, re-download even if cached.

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

        temp_path = f"{filepath}.part"
        if force_download:
            if os.path.exists(filepath):
                os.remove(filepath)
            if os.path.exists(temp_path):
                os.remove(temp_path)

        if not os.path.exists(filepath):
            # Download the file if it doesn't already exist locally
            try:
                existing_bytes = (
                    os.path.getsize(temp_path)
                    if os.path.exists(temp_path)
                    else 0
                )
                headers = (
                    {"Range": f"bytes={existing_bytes}-"}
                    if existing_bytes
                    else {}
                )
                response = requests.get(
                    url, stream=True, headers=headers
                )
                if existing_bytes and response.status_code == 416:
                    response.close()
                    os.remove(temp_path)
                    existing_bytes = 0
                    headers = {}
                    response = requests.get(
                        url, stream=True, headers=headers
                    )
                if existing_bytes and response.status_code != 206:
                    response.close()
                    if os.path.exists(temp_path):
                        os.remove(temp_path)
                    existing_bytes = 0
                    headers = {}
                    response = requests.get(
                        url, stream=True, headers=headers
                    )
                response.raise_for_status()  # Check for HTTP errors
                total_size = int(
                    response.headers.get("Content-Length", 0)
                )
                if existing_bytes and response.status_code == 206:
                    content_range = response.headers.get(
                        "Content-Range"
                    )
                    if content_range and "/" in content_range:
                        total_size = int(
                            content_range.split("/")[-1]
                        )
                    elif total_size:
                        total_size += existing_bytes
                downloaded = existing_bytes
                next_report_percent = 0
                report_step = 5  # percent
                report_bytes_step = 50 * 1024 * 1024
                next_report_bytes = report_bytes_step
                label = (
                    f"{progress_label}: " if progress_label else ""
                )
                if show_progress and not total_size:
                    print(
                        f"{label}starting download (size unknown)"
                    )
                file_mode = "ab" if existing_bytes else "wb"
                with open(temp_path, file_mode) as out_file:
                    for chunk in response.iter_content(
                        chunk_size=1024 * 1024
                    ):
                        if not chunk:
                            continue
                        out_file.write(chunk)
                        downloaded += len(chunk)
                        if show_progress:
                            if total_size:
                                percent = int(
                                    downloaded * 100 / total_size
                                )
                                if percent >= next_report_percent:
                                    print(
                                        f"{label}{percent}% ({downloaded}/{total_size} bytes)"
                                    )
                                    next_report_percent = (
                                        percent + report_step
                                    )
                            elif downloaded >= next_report_bytes:
                                print(
                                    f"{label}{downloaded} bytes"
                                )
                                next_report_bytes += report_bytes_step
                if total_size and downloaded < total_size:
                    raise Exception(
                        "Download incomplete: expected "
                        f"{total_size} bytes, got {downloaded} bytes."
                    )
                os.replace(temp_path, filepath)
            except requests.RequestException as e:
                raise Exception(
                    f"Failed to download the file from {url}. Error: {e}"
                )
            except Exception:
                raise
    else:
        if not os.path.exists(url_or_filepath):
            raise FileNotFoundError(
                f"The file does not exist: {url_or_filepath}"
            )
        filepath = url_or_filepath

    return filepath
