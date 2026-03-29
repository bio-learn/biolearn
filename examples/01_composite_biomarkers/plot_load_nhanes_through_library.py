"""
Loading NHANES Data Through the biolearn Library
=================================================

This example shows the recommended way to load NHANES blood exam data into the
library. NHANES is registered in :class:`~biolearn.data_library.DataLibrary`,
so it loads the same way as every other dataset and returns a
:class:`~biolearn.data_library.GeoData` you can pass straight to a model.
"""

#############################################################################
# Load NHANES 2010 through DataLibrary
# ---------------------------------------
from biolearn.data_library import DataLibrary

data_source = DataLibrary().get("NHANES_2010")
data = data_source.load()

#############################################################################
# Inspect the clinical biomarkers
# ---------------------------------------
# ``data.clinical`` has one row per participant and one column per biomarker.
print("Clinical biomarkers available:")
print(list(data.clinical.columns))

print("\nFirst few rows:")
print(data.clinical.head())

#############################################################################
# Inspect the metadata
# ---------------------------------------
# Age, sex, and mortality fields live on ``data.metadata`` so they stay
# separate from the biomarker measurements.
print("\nMetadata columns:")
print(list(data.metadata.columns))

print("\nFirst few rows of metadata:")
print(data.metadata.head())

#############################################################################
# Plot the age distribution
# ---------------------------------------
import matplotlib.pyplot as plt

data.metadata["age"].hist(bins=30)
plt.xlabel("Age")
plt.ylabel("Number of participants")
plt.title("NHANES 2010: age distribution")
