import yaml
from importlib import resources

def test_sex_values_are_strings():
    lib = yaml.safe_load(resources.files("biolearn.data").joinpath("library.yaml").read_text())
    items = lib.get("datasets", lib.get("items", lib))
    if isinstance(items, list):                       # array layout
        items = {e["id"]: e for e in items}
    valid = {"female", "male", "unknown"}
    for entry in items.values():
        sex = entry.get("metadata", {}).get("sex", entry.get("sex"))
        if sex is None:
            continue
        assert isinstance(sex, str) and sex.lower() in valid
