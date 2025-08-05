"""
Hurdle Inflammage API Example
=============================

This example shows how to use Hurdle Bio's Inflammage model through their API.
"""

import os
import pandas as pd
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import DataLibrary

# Step 1: Set up your API key
# Get your API key from: https://dashboard.sandbox.hurdle.bio/register/partner
# Then either:
# - Set environment variable: export HURDLE_API_KEY="your_key_here"
# - Or pass directly to the model (shown below)

# Step 2: Load some methylation data
print("Loading sample data...")
data = DataLibrary().get("BoAChallengeData").load()

# Use a small subset for this example
sample_data = data.dnam.iloc[:, :5]  # First 5 samples
print(f"Using {sample_data.shape[1]} samples")

# Step 3: Initialize the model
gallery = ModelGallery()

# If you have your API key in environment variable:
# model = gallery.get("HurdleInflammage")

# Or provide it directly:
api_key = input("Enter your Hurdle API key: ").strip()
model = gallery.get("HurdleInflammage")
model.api_key = api_key

# Step 4: Make predictions
# Note: You will be asked for consent before data is sent
print("\nMaking predictions...")
try:
    predictions = model.predict(sample_data)

    print("\nResults:")
    for sample, age in predictions.items():
        print(f"{sample}: {age:.1f} years")

except ValueError as e:
    print(f"Error: {e}")
except Exception as e:
    print(f"API Error: {e}")

# Optional: If you have metadata with chronological ages
if hasattr(data, 'metadata') and 'age' in data.metadata.columns:
    # Compare with chronological age
    chrono_ages = data.metadata.loc[predictions.index, 'age']

    print("\nAge Acceleration:")
    for sample in predictions.index:
        if sample in chrono_ages.index:
            inflamm_age = predictions[sample]
            chrono_age = chrono_ages[sample]
            delta = inflamm_age - chrono_age
            print(f"{sample}: {delta:+.1f} years")
