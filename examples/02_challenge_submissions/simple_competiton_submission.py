"""
Building a competition submission using an existing model
=============================================================

This example shows you how to generate a submission for the community warmup using an existing model
"""

#############################################################################
# Loading up the data for the competition
# ---------------------------------------
import pandas as pd

#Download the data file for the warmup challenge linked here https://www.synapse.org/#!Synapse:syn52966292/wiki/625231
DOWNLOADED_DATA_FILE_PATH="ADD YOUR PATH HERE"
challenge_data = pd.read_csv(DOWNLOADED_DATA_FILE_PATH, index_col=0)
challenge_data

#############################################################################
# Use the Lin model to predict the age
# ------------------------------------------------------
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import GeoData
#Lin scored the best so lets use that
data = GeoData(None, challenge_data)
results = ModelGallery().get("Lin").predict(data)
results

############################################################################################################################
# Save a csv in the correct format for submission
# --------------------------------------------------------------------------------------------------------------------------
#Now we save this to a csv in the expected format. 
submission = results.rename_axis("sampleId")
submission = submission.rename(columns={"Predicted": "predictedAge"})
submission.to_csv("lin_submission.csv")
#You can submit this file to the community warmup challenge by following the instructions below
submission