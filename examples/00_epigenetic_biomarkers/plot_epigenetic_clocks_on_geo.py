"""
"Epigenetic Clocks" in GEO Data
=============================================================

This example loads a DNA Methylation data from GEO, calculates multiple epigenetic clocks, and plots them against chronological age. 
"""

#############################################################################
# First load up some methylation data from GEO using the data library
# ---------------------------------------
from biolearn.data_library import DataLibrary
#Load up GSE41169 blood DNAm data
data_source = DataLibrary().get("GSE41169")
data=data_source.load()

#The data has the methylation data as well as metadata for each subject
methylation_data = data.dnam

#######################################################################################
# Run a simple imputation to replace missing data using the avearges of measured values
# -------------------------------------------------------------------------------------
from biolearn.imputation import impute_from_average
prepared_data = impute_from_average(methylation_data)

######################################################################################
# Now run three different clocks on the dataset to produce epigenetic clock ages
# ------------------------------------------------------------------------------------

from biolearn.clock_gallery import ClockGallery
gallery = ClockGallery()
horvath_results = gallery.get("Horvathv1").predict(prepared_data)
hannum_results = gallery.get("Hannum").predict(prepared_data)
phenoage_results = gallery.get("PhenoAge").predict(prepared_data)



##########################################################################################################
# Finally exctract the age data from the metadata from GEO and plot the results using seaborn
# --------------------------------------------------------------------------------------------------------
import seaborn as sn
import pandas as pd

actual_age = data.metadata['age']
plot_data = pd.DataFrame({
    'Horvath': horvath_results,
    'Hannum': hannum_results,
    'PhenoAge': phenoage_results,
    "Age": actual_age
})
plot_data.index=plot_data['Age']

sn.relplot(data=plot_data, kind="scatter");