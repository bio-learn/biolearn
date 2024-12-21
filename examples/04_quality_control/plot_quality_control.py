"""
Quality Control using GEO Datasets
=============================================================

This example demonstrates the built-in quality control plotting function to visualize 
the distribution of sample deviations from the population mean in a ridge density plot.
"""

################################################################################
# Create a dictionary of datasets id's to datasets
# ------------------------------------------------------------------------------

from biolearn.data_library import DataLibrary, GeoData

library = DataLibrary()
dataset_ids = ["GSE112618","GSE110554","GSE41169","GSE52588"]
datasets = {id:library.get(id).load() for id in dataset_ids}

#############################################################################
# Generates a quality control report for the genomic data containing both 
# detailed methylation data, a summary, and a detailed section for missing percentages per site.
# --------------------------------------------------------------------------

[dataset.quality_report().show() for dataset in datasets.values()]

#############################################################################
# Visualize the distribution of sample deviations from the population mean in a ridge density plot
# --------------------------------------------------------------------------

GeoData.plot_sample_deviations(datasets=datasets)
