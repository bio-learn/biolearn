"""
Aging Clock/Model visualizations using GEO Datasets & ModelAnalyzer Class
=============================================================

This example demonstrates the built-in aging clock/model function to visualize 
different plots on clock(s)/model(s) predictions.
"""

################################################################################
# Import the necessary classes
# ------------------------------------------------------------------------------
from biolearn.data_library import DataLibrary
from biolearn.model_gallery import ModelGallery
from biolearn.visualize import ModelAnalyzer

################################################################################
# Visualize a correlation matrix across aging clocks/models (using plot_clock_correlation_matrix method)
# ------------------------------------------------------------------------------

# Load an appropriate GEO dataset for the models
data = DataLibrary().get("GSE112618").load()

# Instatiate ModelAnalyzer class
model_analyzer = ModelAnalyzer()

modelnames = ["Horvathv1", "Hannum", "PhenoAge", "DunedinPACE", "Lin", "Zhang_10"]
models = [ModelGallery().get(names) for names in modelnames]
model_analyzer.plot_clock_correlation_matrix(
    models=models,
    data=data,
)

################################################################################
# Visualize clock/model chronological age deviations across samples in a heatmap (using plot_clock_deviation_heatmap method)
# ------------------------------------------------------------------------------
modelnames = ["Horvathv1", "Hannum", "PhenoAge", "DunedinPACE", "Lin", "Zhang_10"]
models = [ModelGallery().get(names) for names in modelnames]
model_analyzer.plot_clock_deviation_heatmap(
    models=models,
    data=data,
)

################################################################################
# Visualize aging clock/model predictions against chronological age (using plot_age_predictions method)
# ------------------------------------------------------------------------------
data2 = DataLibrary().get("GSE41169").load()

modelnames = ["Horvathv1", "Hannum", "PhenoAge", "Lin"]
models = [gallery.get(name) for name in modelnames]
model_analyzer.plot_age_predictions(
    models=models,
    data=data2,
)

################################################################################
# Visualize model predictions against it's corresponding health outcome (using plot_clock_deviation_heatmap method)
# ------------------------------------------------------------------------------

# Load an appropriate GEO dataset for the corresponding model
down_syndrome_data = DataLibrary().get("GSE52588").load()

model = [ModelGallery().get("DownSyndrome")]
model_analyzer.plot_health_outcome(
    models=model,
    data=down_syndrome_data,
    # Provide the health outcome column name
    health_outcome_col="disease_state",
)