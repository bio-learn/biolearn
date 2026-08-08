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
    Runs predictions using a collection of models from the model gallery specified in predictors_dict and returns a DataFrame with the results.

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


def calculate_c_index(
    data,
    predictor_results,
    ci_bootstrap_samples=0,
    seed=42,
    adjust_for_age=False,
):
    """
    Calculates the C-index for each predictor in the predictor_results DataFrame.

    Args:
        data (Dataset): A Dataset object containing metadata with columns:
            'dead' - boolean indicating if the subject is dead
            'years_until_death' - time until death or censoring
        predictor_results (pd.DataFrame): DataFrame containing predictor results. Columns are the names of the predictors, and rows are IDs from data.
        ci_bootstrap_samples (int): When > 0, also report a 95% percentile-bootstrap
            confidence interval for each C-index (columns 'CI95_low'/'CI95_high'),
            resampling subjects with replacement this many times. A C-index of 0.75
            from 100 deaths and from 5,000 deaths support very different conclusions;
            the interval makes that difference visible.
        seed (int): Seed for the bootstrap resampling, so reported intervals are
            reproducible.
        adjust_for_age (bool): When True, the C-index is computed on the residuals
            of each predictor after regressing out chronological age (metadata
            column 'age'), mirroring the standardization used in
            calculate_mortality_hazard_ratios. Chronological age alone predicts
            mortality, so an unadjusted C-index rewards a clock merely for
            correlating with age; the adjusted value asks what the clock adds
            beyond it.

    Returns:
        pd.DataFrame: A DataFrame containing C-index values for each predictor,
        plus CI columns when ci_bootstrap_samples > 0.
    """
    # Merge predictor results with metadata
    analysis_df = pd.merge(
        predictor_results, data.metadata, left_index=True, right_index=True
    )

    # Remove rows with missing 'dead' or 'years_until_death' values
    analysis_df = analysis_df.dropna(subset=["dead", "years_until_death"])
    if adjust_for_age:
        analysis_df = analysis_df.dropna(subset=["age"])

    def compute_scores(frame, clock):
        predictor_values = frame[clock].astype(float)
        if adjust_for_age:
            age = frame["age"].astype(float)
            slope, intercept = np.polyfit(age, predictor_values, 1)
            predictor_values = predictor_values - (slope * age + intercept)
        return predictor_values

    c_index_values = []
    ci_lows = []
    ci_highs = []

    for clock in predictor_results.columns:
        predictor_values = compute_scores(analysis_df, clock)

        # Calculate the C-index directly
        c_index = concordance_index(
            event_times=analysis_df["years_until_death"],
            predicted_scores=-predictor_values,  # Negative if higher scores indicate higher risk
            event_observed=analysis_df["dead"],
        )
        c_index_values.append(c_index)

        if ci_bootstrap_samples > 0:
            rng = np.random.default_rng(seed)
            boot_values = []
            n = len(analysis_df)
            for _ in range(ci_bootstrap_samples):
                idx = rng.integers(0, n, size=n)
                sample = analysis_df.iloc[idx]
                if sample["dead"].astype(bool).sum() == 0:
                    continue
                boot_values.append(
                    concordance_index(
                        event_times=sample["years_until_death"],
                        predicted_scores=-compute_scores(sample, clock),
                        event_observed=sample["dead"],
                    )
                )
            ci_lows.append(np.quantile(boot_values, 0.025))
            ci_highs.append(np.quantile(boot_values, 0.975))

    # Create a DataFrame with the results
    results_df = pd.DataFrame(
        {
            "Clock": predictor_results.columns,
            "C_index": c_index_values,
        }
    )
    if ci_bootstrap_samples > 0:
        results_df["CI95_low"] = ci_lows
        results_df["CI95_high"] = ci_highs

    return results_df


def calculate_mortality_hazard_ratios(data, predictor_results):
    """
    Calculates mortality hazard ratios for predictor results using a Cox Proportional Hazards model.

    Args:
        data (GeoData): GeoData object. The metadata must contain the following columns
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
