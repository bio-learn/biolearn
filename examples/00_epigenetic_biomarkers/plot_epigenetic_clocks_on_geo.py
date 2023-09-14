"""
"Epigenetic Clocks" in GEO Data
=============================================================

This example loads a DNA Methylation data from GEO, calculates multiple epigenetic clocks, and plot them against chronological age. 
"""

#############################################################################
# Loading a DNAm GEO data
# ---------------------------------------
from biolearn.data_library import DataLibrary
#Load up GSE41169 blood DNAm data
data_source = DataLibrary().get("GSE41169")
data=data_source.load()

#The data has the methylation data as well as metadata for each subject
df = data.dnam

######################################################################################
# Calculate "biological age" based on Horvath, Hannum, and PhenoAge epigenetic clocks
# -------------------------------------------------------------------------------------
from biolearn.clock import horvathv1, hannum, phenoage
horvath_results = horvathv1(df)
hannum_results = hannum(df)
phenoage_results = phenoage(df)

actual_age = data.metadata['age']

##########################################################################################################
# Plot biological ages against chronological age
# --------------------------------------------------------------------------------------------------------
import seaborn as sn
import pandas as pd
plot_data = pd.DataFrame({
    'Horvath': horvath_results,
    'Hannum': hannum_results,
    'PhenoAge': phenoage_results,
    "Age": actual_age
})
plot_data.index=plot_data['Age']

sn.relplot(data=plot_data, kind="scatter");