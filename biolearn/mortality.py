import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter

from biolearn.model_gallery import ModelGallery


def run_predictions(data, predictors_dict):
    """
    Runs predictions using a collection of models specified in predictors_dict and returns a DataFrame with the results.

    Args:
        data (GeoData): GeoData object used for predictions
        predictors_dict (dict): A dictionary where keys are model names and values are lists of column names
                                in the prediction output to be used. ex ("Horvathv1": "Predicted")

    Returns:
        pd.DataFrame: A DataFrame containing predictions from each model.
    """
    # DataFrame to store results
    results_df = pd.DataFrame()
    gallery = ModelGallery()

    # Loop through each model and make predictions
    for model_name, keys in predictors_dict.items():
        model = gallery.get(model_name)
        prediction = model.predict(data)
        results_df[model_name] = prediction[keys]

    return results_df


def calculate_mortality_hazard_ratios(data, predictor_results):
    """
    Calculates mortality hazard ratios for predictor results using a Cox Proportional Hazards model.

    Args:
        data (GeoData): GeoData object. The metadata must contain the followin columns
             'age' - age in years
             'dead' - 0 for alive, 1 for dead
             'years_until_death' - if dead this should be years between sample collection and death. Otherwise years between sample collection and last known contact with live subject
        predictor_results (pd.DataFrame): The DataFrame containing predictor results. Columns are the name of the predictor model and rows must be ids from data

    Returns:
        pd.DataFrame: A DataFrame containing hazard ratios, confidence intervals, and p-values for each predictor.
    """
    analysis_df = pd.merge(
        predictor_results, data.metadata, left_index=True, right_index=True
    )

    # Standardize the clock values
    for clock in predictor_results.columns:
        analysis_df[clock] = (
            analysis_df[clock] - analysis_df[clock].mean()
        ) / analysis_df[clock].std()

    # Remove rows where 'dead' column is null
    analysis_df = analysis_df.dropna(subset=["dead"])

    hazard_ratios = []
    ci_lower_list = []
    ci_upper_list = []
    p_values = []

    for clock in predictor_results.columns:
        cph = CoxPHFitter()
        cph.fit(
            analysis_df[["age", clock, "years_until_death", "dead"]],
            duration_col="years_until_death",
            event_col="dead",
        )
        hazard_ratios.append(cph.hazard_ratios_[clock])
        ci_lower, ci_upper = cph.confidence_intervals_.loc[clock]
        ci_lower_list.append(np.exp(ci_lower))
        ci_upper_list.append(np.exp(ci_upper))
        p_value = cph.summary.loc[clock, "p"]
        p_values.append(p_value)

    results_df = pd.DataFrame(
        {
            "Clock": predictor_results.columns,
            "HR": hazard_ratios,
            "CI_lower": ci_lower_list,
            "CI_upper": ci_upper_list,
            "P_value": p_values,
        }
    )

    return results_df


def plot_hazard_ratios(hazard_ratio_data):
    """
    Plots hazard ratios from the provided data in a forest plot.

    Args:
        hazard_ratio_data (pd.DataFrame): A DataFrame containing hazard ratios, confidence intervals, and p-values
                                          for each predictor.

    Returns:
        None
    """
    sorted_data = hazard_ratio_data.sort_values("HR", ascending=True)

    # Generate colors from the "rocket" palette
    num_clocks = len(sorted_data)
    colors = sns.color_palette("rocket", num_clocks)[::-1]

    # Create the forest plot with the updated colors
    fig, ax = plt.subplots(figsize=(16, 6))

    for i, (index, row) in enumerate(sorted_data.iterrows()):
        ax.plot(
            [row["CI_lower"], row["CI_upper"]],
            [i, i],
            color=colors[i],
            linewidth=2,
            alpha=0.7,
        )
        ax.scatter(row["HR"], i, color=colors[i], s=100, zorder=3)
        ax.annotate(
            f"HR = {row['HR']:.2f} ({row['CI_lower']:.2f}-{row['CI_upper']:.2f}), P = {row['P_value']:.2e}",
            xy=(row["CI_upper"], i),
            xytext=(6, 0),
            textcoords="offset points",
            va="center",
            color=colors[i],
        )

    ax.set_yticks(range(len(sorted_data)))
    ax.set_yticklabels(sorted_data["Clock"])
    ax.set_xlabel("Hazard Ratio (95% CI) per SD increase, adjusted for age")
    ax.set_title("Mortality Hazard Ratio From Clock Predictions")
    ax.axvline(x=1, color="black", linestyle="--", linewidth=1)

    plt.tight_layout()
    plt.show()
