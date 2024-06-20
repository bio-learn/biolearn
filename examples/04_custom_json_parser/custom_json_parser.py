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


library = DataLibrary(get_data_file("new_all_geo_ids.yaml"))
print(library.cache.path)

data_source = library.get("GSE56046")
data = data_source.load()

print("\nMetadata:")
print("-" * 40)
print(data.metadata)

print("\nMethylation data:")
print("-" * 40)
print(data.dnam)

if (data.dnam is None or data.dnam.shape[0] == 0):
    print("\nList Supplementary Files:")
    print("-" * 40)
    if "supplementary_files" in data_source.data:
        for suppl_file in data_source.data['supplementary_files']:
            print(suppl_file['name'], suppl_file['size'], suppl_file['links'])
