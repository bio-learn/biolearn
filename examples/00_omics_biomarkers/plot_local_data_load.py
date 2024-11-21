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

file_path = get_test_data_file("example_data.csv")
data = GeoData.from_methylation_matrix(file_path)

data.dnam


######################################################################################
# You can now use it like any other Biolearn dataset
# ------------------------------------------------------------------------------------

data.quality_report().show()