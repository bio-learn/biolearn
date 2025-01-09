"""
Quality control visualization using GEO datasets
=============================================================

This example demonstrates the built-in quality control plotting function to visualize 
the distribution of sample deviations from the population mean in a ridge density plot.
"""

################################################################################
# Import required classes and functions
# ------------------------------------------------------------------------------
from biolearn.data_library import DataLibrary, GeoData
from biolearn.visualize import plot_sample_deviations

################################################################################
# Create a dictionary of dataset display names to GeoData objects
# ------------------------------------------------------------------------------
library = DataLibrary()
dataset_ids = ["GSE112618", "GSE110554", "GSE41169", "GSE52588"]
datasets = {id: library.get(id).load() for id in dataset_ids}

################################################################################
# Generate a quality control report for each dataset
# ------------------------------------------------------------------------------
[dataset.quality_report().show() for dataset in datasets.values()]

################################################################################
# Visualize the distribution of sample deviations from the population mean
# ------------------------------------------------------------------------------

# Use the `plot_sample_deviations` function to generate a ridge density plot
plot_sample_deviations(datasets=datasets)
