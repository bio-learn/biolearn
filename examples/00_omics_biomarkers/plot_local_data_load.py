"""
Local Data Loading
=============================================================

This example loads data from a local file
"""

#############################################################################
# Load up a local data file
# ---------------------------------------------------------------------------
from biolearn.data_library import GeoData
from biolearn.util import get_test_data_file

#Files formatted as described in the standard https://bio-learn.github.io/methylation-standard.html will load correctly. 
#Load will search for files names [name]_metadata.csv and [name]_methylation_part0.csv

file_path = get_test_data_file("")
data = GeoData.load_csv(file_path, "example")

data.dnam

#############################################################################
# Metadata is loaded if available
# ---------------------------------------------------------------------------
data.metadata

######################################################################################
# You can now use it like any other Biolearn dataset
# ------------------------------------------------------------------------------------

data.quality_report().show()