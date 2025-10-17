"""
Hurdle InflammAge API Example
==============================

This example demonstrates using Hurdle Bio's InflammAge model through their API.
InflammAge is an inflammation-based biological age clock.

Before running this example:
1. Register at https://dashboard.hurdle.bio/register
2. Navigate to Developer > API keys
3. Generate and save your API key
4. Set it as an environment variable: export HURDLE_API_KEY="your_key_here"
"""

import os
import pandas as pd
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import DataLibrary

print("Loading sample data...")
data = DataLibrary().get("BoAChallengeData").load()

# Use first 5 samples for quick testing
sample_data = data.dnam.iloc[:, :5]
print(f"Using {sample_data.shape[1]} samples with {sample_data.shape[0]} CpG sites")

gallery = ModelGallery()

# Check if API key is set
if not os.environ.get("HURDLE_API_KEY"):
    print("\nNo HURDLE_API_KEY environment variable found.")
    api_key = input("Enter your Hurdle API key: ").strip()
    if api_key:
        os.environ["HURDLE_API_KEY"] = api_key
    else:
        print("No API key provided. Exiting.")
        exit(1)

# Load the model
model = gallery.get("HurdleInflammAge")

print("\nMaking predictions...")
print("Note: You will be prompted to consent to sending data to Hurdle's servers.")

try:
    predictions = model.predict(sample_data)

    print("\nInflammAge Results:")
    print("-" * 40)
    for sample, age in predictions.items():
        print(f"{sample}: {age:.1f} years")

    # If metadata is available, show age acceleration
    if hasattr(data, "metadata") and "age" in data.metadata.columns:
        print("\nAge Acceleration (InflammAge - Chronological Age):")
        print("-" * 40)
        for sample in predictions.index:
            if sample in data.metadata.index:
                inflamm_age = predictions[sample]
                chrono_age = data.metadata.loc[sample, "age"]
                delta = inflamm_age - chrono_age
                print(f"{sample}: {delta:+.1f} years")

except ValueError as e:
    print(f"Error: {e}")
except Exception as e:
    print(f"API Error: {e}")
    print("\nTroubleshooting:")
    print("- Verify your API key is correct")
    print("- Check your internet connection")
    print("- Ensure you have access to https://api.hurdle.bio")
