import os
from biolearn.data_library import DataLibrary


def get_data_file(relative_path):
    script_dir = os.path.dirname(
        __file__
    )  # get the directory of the current script
    data_file_path = os.path.join(
        script_dir, "data", relative_path
    )  # build the path to the data file
    return data_file_path


data = DataLibrary(get_data_file("new_all_geo_ids.yaml")).get("GSE28094").load()
print(data.metadata)
print("-" * 40)
print(data.dnam)
