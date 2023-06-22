"""
Plotting survival curves with NHANES data
=============================================================

This example loads up data from NHANES and demonstrates plotting survival curves 
including with a PhenoAge clock
"""

#############################################################################
# Needed imports
# ---------------------------------------
import pandas as pd
import numpy as np
from warnings import simplefilter

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

from biolearn.load import load_nhanes

#############################################################################
# Load up data from NHANES
# ---------------------------------------
year = 2010
df = load_nhanes(year)
df["LBXCRP"] = np.log(df["LBXCRP"])
df["years_until_death"] = df["months_until_death"] / 12

#############################################################################
# Basic survival plot
# ---------------------------------------
T = df.years_until_death
E = df.is_dead
kmf = KaplanMeierFitter()
kmf.fit(T, E)
kmf.plot()
plt.ylabel("Survival")
plt.xlabel("Years")

#############################################################################
# Plot survival by sex
# ---------------------------------------
ax = plt.subplot()
groups = df["sex"]
ix = groups == 2
kmf.fit(T[ix], E[ix], label="Female")
ax = kmf.plot_survival_function(ax=ax)
kmf.fit(T[~ix], E[~ix], label="Male")
ax = kmf.plot_survival_function()
plt.ylabel("Survival")
plt.xlabel("Years")

#############################################################################
# Plot survival by glucose level
# ---------------------------------------
ax = plt.subplot()
groups = df["glucose"]
ix = groups < 5.5
kmf.fit(T[ix], E[ix], label="glucose")
ax = kmf.plot_survival_function(ax=ax)
kmf.fit(T[~ix], E[~ix], label="glucose")
ax = kmf.plot_survival_function()
plt.ylabel("Survival")
plt.xlabel("Years")

#############################################################################
# Calculate "biological age" using the phenoage clock
# ------------------------------------------------------
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940111/
pheno_coefs = {
    "age": 0.0804,
    "LBDSALSI": -0.034,
    "LBDSCRSI": 0.0095,
    "glucose": 0.1953,
    "LBXCRP": 0.0954,
    "LBXLYPCT": -0.012,
    "LBXMCVSI": 0.0268,
    "LBXRDW": 0.3356,
    "LBXSAPSI": 0.00188,
    "LBXWBCSI": 0.0554,
}

constant = -19.9067
gamma = 0.0077
cs = [141.50225, -0.00553, 0.090165]
df["pheno"] = df.apply(
    lambda x: sum([x[c] * pheno_coefs[c] for c in pheno_coefs.keys()]), axis=1
)
df["mortality_score"] = 1 - np.exp(
    (-np.exp(constant + df["pheno"]) * ((np.exp(120 * gamma) - 1) / gamma))
)
df["phenotypic_age"] = (
    cs[0] + np.log(cs[1] * np.log(1 - df["mortality_score"])) / cs[2]
)

#############################################################################
# Show relation between biological and chronological age
# ---------------------------------------------------------------
df.plot.scatter("age", "phenotypic_age")

############################################################################################
# Plot survival curve for people with a biological age younger vs older than chronological
# -----------------------------------------------------------------------------------------
df["biologically_older"] = df["phenotypic_age"] > df["age"]
ax = plt.subplot()
groups = df["biologically_older"]
ix = groups == 0
kmf.fit(T[ix], E[ix], label="Biologically younger")
ax = kmf.plot_survival_function(ax=ax)
kmf.fit(T[~ix], E[~ix], label="Biologically older")
ax = kmf.plot_survival_function()
plt.ylabel("Survival")
plt.xlabel("Years")
