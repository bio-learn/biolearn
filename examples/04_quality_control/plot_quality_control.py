"""
Quality Control using GEO Datasets
=============================================================

This example demonstrates the built-in quality control plotting function to visualize 
the distribution of sample deviations from the population mean in a ridge density plot.
"""

################################################################################
# Load GEO methylation dataset(s) to use as a reference for quality control
# ------------------------------------------------------------------------------

from biolearn.data_library import DataLibrary, GeoData

# Load the reference dataset and compute averages for imputation
dataset1 = DataLibrary().get('GSE112618').load()
dataset2 = DataLibrary().get('GSE110554').load()
dataset3 = DataLibrary().get('GSE41169').load()
dataset4 = DataLibrary().get('GSE52588').load()

#############################################################################
# Put the datasets and dataset ids in separate lists
# --------------------------------------------------------------------------

datasets = [dataset1,dataset2,dataset3,dataset4]
dataset_ids = ["GSE112618","GSE110554","GSE41169","GSE52588"]

#############################################################################
# Generates a quality control report for the genomic data containing both 
# detailed methylation data, a summary, and a detailed section for missing percentages per site.
# --------------------------------------------------------------------------

dataset1.quality_report().show()
dataset2.quality_report().show()
dataset3.quality_report().show()
dataset4.quality_report().show()

#############################################################################
# Visualize the distribution of sample deviations from the population mean in a ridge density plot
# --------------------------------------------------------------------------

GeoData.plot_sample_deviations(datasets=datasets, dataset_ids=dataset_ids)
