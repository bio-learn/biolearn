"""
"Epigenetic Clocks" in GEO Data
=============================================================

This example loads a DNA Methylation data from GEO, calculates multiple epigenetic clocks, and plot them against chronological age. 
"""

#############################################################################
# Loading a DNAm GEO data
# ---------------------------------------
from biolearn.load import load_dnam
#GSE41169 blood DNAm data
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41169/matrix/GSE41169_series_matrix.txt.gz'
df=load_dnam(dnam_file=url,id_row=32,age_row=46,skiprows=72)

#############################################################################
# Calculate "biological age" based on Horvath, Hannum, and PhenoAge epigenetic clocks
# ------------------------------------------------------
from biolearn.clock import horvath_clock, hannum_clock, phenoage_clock
df['horvath'] = horvath_clock(df)
df['hannum'] = hannum_clock(df)
df['phenoage'] = phenoage_clock(df)

##########################################################################################################
# Plot biological ages against chronological age
# --------------------------------------------------------------------------------------------------------
import seaborn as sn
df.index=df['age']
sn.relplot(data=df[['horvath','hannum','phenoage']], kind="scatter");