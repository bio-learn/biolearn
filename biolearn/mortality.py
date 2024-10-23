import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index

from biolearn.model_gallery import ModelGallery
import warnings


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
        try:
            model = gallery.get(model_name)
            prediction = model.predict(data)
            results_df[model_name] = prediction[keys]
        except Exception as e:
            # Catch any errors, issue a warning, and continue with the next model
            warnings.warn(
                f"Error running model '{model_name}': {str(e)}", RuntimeWarning
            )
            continue

    return results_df


def calculate_c_index(data, predictor_results):
    """
    Calculates the C-index for each predictor in the predictor_results DataFrame without adjusting for age.

    Args:
        data (Dataset): A Dataset object containing metadata with columns:
            'dead' - boolean indicating if the subject is dead
            'years_until_death' - time until death or censoring
        predictor_results (pd.DataFrame): DataFrame containing predictor results. Columns are the names of the predictors, and rows are IDs from data.

    Returns:
        pd.DataFrame: A DataFrame containing C-index values for each predictor.
    """
    # Merge predictor results with metadata
    analysis_df = pd.merge(
        predictor_results, data.metadata, left_index=True, right_index=True
    )

    # Remove rows with missing 'dead' or 'years_until_death' values
    analysis_df = analysis_df.dropna(subset=["dead", "years_until_death"])

    c_index_values = []

    for clock in predictor_results.columns:
        # Ensure predictor values are numeric
        predictor_values = analysis_df[clock].astype(float)

        # Calculate the C-index directly
        c_index = concordance_index(
            event_times=analysis_df["years_until_death"],
            predicted_scores=-predictor_values,  # Negative if higher scores indicate higher risk
            event_observed=analysis_df["dead"],
        )
        c_index_values.append(c_index)

    # Create a DataFrame with the results
    results_df = pd.DataFrame(
        {
            "Clock": predictor_results.columns,
            "C_index": c_index_values,
        }
    )

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


def calculate_age_adjusted_c_index(data, predictor_results):
    """
    Calculates the C-index for each predictor in the predictor_results DataFrame, adjusted for age.

    Args:
        data (Dataset): A Dataset object containing metadata with columns:
            'dead' - boolean indicating if the subject is dead
            'years_until_death' - if dead this should be years between sample collection and death. Otherwise years between sample collection and last known contact with live subject
            'age' - age of the subject at sample collection
        predictor_results (pd.DataFrame): The DataFrame containing predictor results. Columns are the name of the predictor model and rows must be ids from data

    Returns:
        pd.DataFrame: A DataFrame containing C-index values for each predictor.
    """
    analysis_df = pd.merge(
        predictor_results, data.metadata, left_index=True, right_index=True
    )

    # Remove rows where 'dead' column is null
    analysis_df = analysis_df.dropna(subset=["dead", "age"])

    c_index_values = []

    for clock in predictor_results.columns:
        cph = CoxPHFitter()
        cph.fit(
            analysis_df[[clock, "age", "years_until_death", "dead"]],
            duration_col="years_until_death",
            event_col="dead",
        )
        c_index = concordance_index(
            analysis_df["years_until_death"],
            -cph.predict_partial_hazard(analysis_df),
            analysis_df["dead"],
        )
        c_index_values.append(c_index)

    results_df = pd.DataFrame(
        {
            "Clock": predictor_results.columns,
            "C_index": c_index_values,
        }
    )

    return results_df


def calculate_log_rank_test(data, predictor_results):
    """
    Calculates the log-rank test for each predictor in the predictor_results DataFrame, adjusted for age.

    Args:
        data (Dataset): A Dataset object containing metadata with columns:
            'dead' - boolean indicating if the subject is dead
            'years_until_death' - if dead this should be years between sample collection and death. Otherwise years between sample collection and last known contact with live subject
            'age' - age of the subject at sample collection
        predictor_results (pd.DataFrame): The DataFrame containing predictor results. Columns are the name of the predictor model and rows must be ids from data

    Returns:
        pd.DataFrame: A DataFrame containing log-rank test statistics and p-values for each predictor.
    """
    analysis_df = pd.merge(
        predictor_results, data.metadata, left_index=True, right_index=True
    )

    # Remove rows where 'dead' column is null
    analysis_df = analysis_df.dropna(subset=["dead", "age"])

    log_rank_stats = []
    p_values = []

    for clock in predictor_results.columns:
        # Fit Cox model with age and clock
        cph = CoxPHFitter()
        cph.fit(
            analysis_df[[clock, "age", "years_until_death", "dead"]],
            duration_col="years_until_death",
            event_col="dead",
        )

        # Calculate age-adjusted risk scores
        risk_scores = cph.predict_partial_hazard(analysis_df)

        # Divide the age-adjusted risk scores into two groups based on median
        median_value = risk_scores.median()
        high_risk = risk_scores > median_value

        # Perform log-rank test
        results = logrank_test(
            analysis_df["years_until_death"][~high_risk],
            analysis_df["years_until_death"][high_risk],
            analysis_df["dead"][~high_risk],
            analysis_df["dead"][high_risk],
        )

        log_rank_stats.append(results.test_statistic)
        p_values.append(results.p_value)

    results_df = pd.DataFrame(
        {
            "Clock": predictor_results.columns,
            "Chi_square": log_rank_stats,
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
