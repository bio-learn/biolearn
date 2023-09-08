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
df['horvath'] = horvathv1(df)
df['hannum'] = hannum(df)
df['phenoage'] = phenoage(df)

##########################################################################################################
# Plot biological ages against chronological age
# --------------------------------------------------------------------------------------------------------
import seaborn as sn
df.index=df['age']
sn.relplot(data=df[['horvath','hannum','phenoage']], kind="scatter");