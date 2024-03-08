"""
Training an ElasticNet model
=============================================================

This example shows you how to train a simple elasticnet model and use it to submit to the challenge
"""

#############################################################################
# Loading up the data for the competition
# ---------------------------------------
from biolearn.data_library import GeoData

#Download the data file for the warmup challenge linked here https://www.synapse.org/#!Synapse:syn52966292/wiki/625231
DOWNLOADED_DATA_FILE_PATH="ADD YOUR PATH HERE"
challenge_data = GeoData.from_methylation_matrix(DOWNLOADED_DATA_FILE_PATH)
challenge_data


#############################################################################
# Load up some training data
# ------------------------------------------------------
from biolearn.data_library import DataLibrary
data = DataLibrary().get("GSE40279").load()
data.metadata


############################################################################################################################
# Narrow down what sites are correlated with age
# --------------------------------------------------------------------------------------------------------------------------
#NOTE: This takes a long time to run
import numpy as np
from sklearn.linear_model import LinearRegression

# Extract data from your 'data' object
X = data.dnam.transpose().values  # Transpose to have samples as rows and cpg sites as columns
y = data.metadata['age'].values

# Parameters for bootstrap and feature selection
n_bootstrap = 20
threshold = 0.05 

# Store count of times each CpG site is deemed significant
cpg_counts = np.zeros(X.shape[1])

# Begin bootstrap iterations
for _ in range(n_bootstrap):
    # Sample with replacement from X, y
    sample_idx = np.random.choice(range(X.shape[0]), size=X.shape[0], replace=True)
    X_sample = X[sample_idx]
    y_sample = y[sample_idx]
    
    # Train model
    model = LinearRegression()
    model.fit(X_sample, y_sample)
    
    # Identify significant CpG sites (based on magnitude of coefficients)
    significant_cpgs = np.where(np.abs(model.coef_) > threshold)[0]
    cpg_counts[significant_cpgs] += 1

# Determine stable CpG sites
stable_cpg_sites = np.where(cpg_counts > n_bootstrap * 0.6)[0]
stable_cpg_names = data.dnam.index[stable_cpg_sites].tolist()

print(f"Stable CpG sites (associated with age in more than 60% of bootstrap samples): {stable_cpg_sites}")

############################################################################################################################
# Seperate data into training and test sets
# --------------------------------------------------------------------------------------------------------------------------
from sklearn.model_selection import train_test_split

df = data.dnam.transpose()
df['age'] = data.metadata['age']
top_sites_df = df[stable_cpg_names]

X = top_sites_df
y = df['age']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

############################################################################################################################
# Train a model using elastic net
# --------------------------------------------------------------------------------------------------------------------------
# Define the model
# alpha is the regularization strength, and l1_ratio defines the mix between L1 and L2
# l1_ratio = 1 is Lasso; l1_ratio = 0 is Ridge.
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error

model = ElasticNet(alpha=0.01, l1_ratio=0.3, max_iter=10000)
# Train the model
model.fit(X_train, y_train)

# Predict and evaluate
y_pred = model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)

print(f"Mean Squared Error on Test Data: {mse}")

############################################################################################################################
# Plot the results to see how good our model is
# --------------------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt

y_pred = model.predict(X_test)

plt.figure(figsize=(10, 8))
plt.scatter(y_test, y_pred, alpha=0.7, edgecolors='w', linewidth=0.5)
plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'k--', lw=3)  # y=x line for reference
plt.xlabel('Actual Age')
plt.ylabel('Predicted Age')
plt.title('Actual Age vs. Predicted Age')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.show()

# Calculate Mean Squared Error (MSE)
mse = np.mean((y_test - y_pred) ** 2)
print(f"Mean Squared Error (MSE): {mse:.4f}")

# Calculate Mean Absolute Error (MAE)
mae = np.mean(np.abs(y_test - y_pred))
print(f"Mean Absolute Error (MAE): {mae:.4f}")

############################################################################################################################
# Run the challenge data through the model
# --------------------------------------------------------------------------------------------------------------------------

pruned_data = challenge_data.dnam.T[stable_cpg_names]
pruned_data = pruned_data.fillna(0)
challenge_results = model.predict(pruned_data)

############################################################################################################################
# Save the results as an output file for submission
# --------------------------------------------------------------------------------------------------------------------------
import pandas as pd

predicted_age_df = pd.DataFrame({
    'predictedAge': challenge_results
}, index=challenge_data.dnam.columns)
predicted_age_df.index.name = 'sampleId'
predicted_age_df