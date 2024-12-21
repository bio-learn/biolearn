"""
DNA Methylation visualizations using the MethylationAnalyzer Class
=============================================================

This example demonstrates the built-in DNA methylation plotting function to visualize 
the DNA methylation against age across datasets in a violin and linear and polynomial regression plots.
"""

################################################################################
# Import the necessary classes
# ------------------------------------------------------------------------------
from biolearn.data_library import DataLibrary, GeoData, 
from biolearn.visualize import MethylationAnalyzer

################################################################################
# Load GEO methylation dataset(s) and set up objects for visualization function use
# ------------------------------------------------------------------------------
# List dataset IDs
dataset_ids = ['GSE112618', 'GSE110554', 'GSE41169', 'GSE52588']

# Load datasets and combine DNAm data and metadata into a list of tuples
datasets = [(data.dnam, data.metadata) for data in [DataLibrary().get(id).load() for id in dataset_ids]]

# Initialize the MethylationAnalyzer class 
# Takes in a list of tuples containiing dnam and metadata, and dataset ids
dnam_analysis = MethylationAnalyzer(datasets, dataset_ids)

# Identify stable CpG sites with low variance
stable_sites = dnam_analysis.identify_stable_cpg_sites(threshold=0.1)
print("Stable CpG Sites:", stable_sites)

#############################################################################
# Compute statistics by age group for the specified CpG site
# --------------------------------------------------------------------------
# Initialize an appropriate CpG site
cpg_site = "cg00000029"

# Compute statistics by age group for the specified CpG site
stats_df, combined_data = dnam_analysis.compute_age_group_stats(cpg_site, age_bins=5)

# Compute regression statistics to plot methylation vs. age
combined_data, statistics = dnam_analysis.compute_methylation_stats(cpg_site)


#############################################################################
# Visualize DNA methylation against age in a violin and linear and polynomial regression plots
# --------------------------------------------------------------------------
# Plot DNA methylation levels by age group and sex across datasets in a violin plot
dnam_analysis.plot_methylation_by_age_sex(combined_data, cpg_site)

# Plot DNA methylation against age across datasets with linear and polynomial regression lines
dnam_analysis.plot_methylation_vs_age(combined_data, statistics, cpg_site)