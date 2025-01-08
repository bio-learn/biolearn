from biolearn.data_library import DataLibrary, GeoData
from biolearn.model_gallery import ModelGallery
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import ttest_ind, pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import r2_score


def identify_stable_cpg_sites(datasets, threshold=0.01):
    """
    Identifies stable CpG sites in each dataset based on variance.
    Parameters:
    - datasets: Dictionary of dataset names to GeoData objects.
    - threshold: Variance threshold to determine stability (default: 0.01).
    Returns:
    - stable_sites_summary: List of dictionaries with stable CpG sites per dataset.
    """
    stable_sites_summary = []
    for dataset_name, geo_data in datasets.items():
        # Compute variance for each CpG site (rows of DNA methylation matrix)
        variances = geo_data.dnam.var(axis=1)
        # Identify CpG sites with variance below the threshold
        stable_sites = variances[variances < threshold].index.tolist()
        # Append results to the summary list
        stable_sites_summary.append(
            {"Dataset": dataset_name, "Stable CpG Sites": stable_sites}
        )
    return stable_sites_summary


def compute_age_group_stats(datasets, cpg_site, age_bins=5):
    """
    Computes statistics (mean, median, std, etc.) for methylation levels by age group and sex.

    Parameters:
    - datasets: Dictionary of dataset names to GeoData objects.
    - cpg_site: CpG site ID to compute statistics for.
    - age_bins: Number of age bins for grouping.

    Returns:
    - stats_df: DataFrame containing computed statistics.
    - combined_data: DataFrame with added 'age_group' column.
    """
    combined_data = []

    # Combine data from all datasets
    for dataset_name, geo_data in datasets.items():
        # Check if the CpG site exists in the dataset
        if cpg_site not in geo_data.dnam.index:
            print(
                f"Warning: CpG site '{cpg_site}' not found in dataset {dataset_name}. Skipping this dataset."
            )
            continue

        # Extract methylation data for the specified CpG site
        methylation_data = geo_data.dnam.loc[cpg_site]
        # Find valid samples with both metadata and methylation data
        valid_samples = geo_data.metadata.index.intersection(
            methylation_data.index
        )
        # Prepare plot data by merging metadata with methylation values
        plot_data = geo_data.metadata.loc[valid_samples, ["age", "sex"]].copy()
        plot_data["methylation"] = methylation_data[valid_samples]
        plot_data["dataset"] = dataset_name
        combined_data.append(plot_data)

    # Concatenate all datasets into a single DataFrame
    combined_data = pd.concat(combined_data).dropna(
        subset=["age", "methylation"]
    )

    # Define bins for age groups with whole numbers
    min_age = int(np.floor(combined_data["age"].min()))
    max_age = int(np.ceil(combined_data["age"].max()))
    bins = np.linspace(min_age, max_age, age_bins + 1).astype(int)
    # Assign age groups based on defined bins
    combined_data["age_group"] = pd.cut(
        combined_data["age"],
        bins=bins,
        include_lowest=True,
        labels=[f"{bins[i]}-{bins[i+1]-1}" for i in range(len(bins) - 1)],
    )

    stats_summary = []

    # Group by age_group to compute statistics for each age range
    for age_group in combined_data["age_group"].unique():
        group_data = combined_data[combined_data["age_group"] == age_group]

        # Separate data by sex within the age group
        male_data = group_data[group_data["sex"] == "2"]["methylation"]
        female_data = group_data[group_data["sex"] == "1"]["methylation"]

        # Calculate statistics
        stats_row = {
            "age_group": age_group,
            "male_mean": male_data.mean(),
            "female_mean": female_data.mean(),
            "male_median": male_data.median(),
            "female_median": female_data.median(),
            "male_std": male_data.std(),
            "female_std": female_data.std(),
            "male_count": male_data.count(),
            "female_count": female_data.count(),
        }

        # Append statistics for the current age group
        stats_summary.append(stats_row)

    # Convert the summary list to a DataFrame for easy viewing
    stats_df = pd.DataFrame(stats_summary)
    return stats_df, combined_data


def plot_methylation_by_age_sex(combined_data, cpg_site):
    """
    Visualizes methylation levels at a specific CpG site by age group and sex using combined data.

    Parameters:
    - combined_data: DataFrame with combined data from all datasets, including 'age_group' column.
    - cpg_site: CpG site ID to visualize.
    """
    plt.figure(figsize=(14, 8))
    # Create a violin plot for methylation levels grouped by age group and sex
    sns.violinplot(
        x="age_group",
        y="methylation",
        hue="sex",
        data=combined_data,
        split=True,
        palette="Set2",
    )
    plt.xlabel("Age Group")
    plt.ylabel(f"Methylation Level at {cpg_site}")
    plt.legend(title="Sex")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


def compute_methylation_stats(datasets, cpg_site, degree=2):
    """
    Computes regression statistics (linear and polynomial) for DNA methylation vs. age.

    Parameters:
    - datasets: Dictionary of dataset names to GeoData objects.
    - cpg_site: CpG site ID to compute statistics for.
    - degree: Degree of the polynomial for the polynomial fit (default: 2).

    Returns:
    - combined_data: DataFrame with combined methylation and age data.
    - statistics: Dictionary containing linear and polynomial regression results.
    """
    combined_data = []

    # Combine data from all datasets
    for dataset_name, geo_data in datasets.items():
        if cpg_site not in geo_data.dnam.index:
            print(
                f"Warning: CpG site '{cpg_site}' not found in dataset {dataset_name}. Skipping this dataset."
            )
            continue

        # Extract methylation data for the CpG site
        methylation_data = geo_data.dnam.loc[cpg_site]
        valid_samples = geo_data.metadata.index.intersection(
            methylation_data.index
        )
        plot_data = geo_data.metadata.loc[valid_samples, ["age"]].copy()
        plot_data["methylation"] = methylation_data[valid_samples]
        plot_data["dataset"] = dataset_name
        combined_data.append(plot_data)

    # Combine data across datasets
    combined_data = pd.concat(combined_data).dropna(
        subset=["age", "methylation"]
    )

    # Prepare data for regression
    ages = combined_data["age"].values.reshape(-1, 1)
    methylation_levels = combined_data["methylation"].values

    # Perform linear regression
    linear_model = LinearRegression()
    linear_model.fit(ages, methylation_levels)
    linear_pred = linear_model.predict(ages)
    r_squared_linear = r2_score(methylation_levels, linear_pred)
    corr, p_value = pearsonr(ages.flatten(), methylation_levels)

    # Perform polynomial regression
    poly = PolynomialFeatures(degree)
    ages_poly = poly.fit_transform(ages)
    poly_model = LinearRegression()
    poly_model.fit(ages_poly, methylation_levels)
    poly_pred = poly_model.predict(ages_poly)
    r_squared_poly = r2_score(methylation_levels, poly_pred)

    # Compile statistics
    statistics = {
        "linear_pred": linear_pred,
        "r_squared_linear": r_squared_linear,
        "corr": corr,
        "p_value": p_value,
        "poly_pred": poly_pred,
        "r_squared_poly": r_squared_poly,
        "degree": degree,
    }

    # Print statistics summary
    print(f"Computed statistics for CpG site: {cpg_site}")
    print(f"Linear Regression R^2: {r_squared_linear:.4f}")
    print(f"Pearson Correlation (r): {corr:.4f}")
    print(f"Pearson Correlation p-value: {p_value:.4e}")
    print(f"Polynomial Regression (Degree {degree}) R^2: {r_squared_poly:.4f}")

    return combined_data, statistics


def plot_methylation_vs_age(combined_data, statistics, cpg_site):
    """
    Plots the DNA methylation level of a specific CpG site against age,
    with both a linear and polynomial regression line.

    Parameters:
    - combined_data: DataFrame with combined methylation and age data.
    - statistics: Dictionary containing regression statistics.
    - cpg_site: CpG site ID to visualize.
    """
    # Extract age and methylation data from the combined dataset
    ages = combined_data["age"].values
    methylation_levels = combined_data["methylation"].values

    plt.figure(figsize=(12, 6))

    # Scatter plot of methylation levels, color-coded by dataset
    sns.scatterplot(
        x=ages,
        y=methylation_levels,
        hue=combined_data["dataset"],
        palette="tab10",
        s=10,
        alpha=0.6,
        legend="full",
    )

    # Linear fit line
    sns.lineplot(
        x=ages, y=statistics["linear_pred"], color="red", linestyle="--"
    )

    # Polynomial fit line
    sns.lineplot(x=ages, y=statistics["poly_pred"], color="green")

    plt.xlabel("Age")
    plt.ylabel("DNA Methylation")
    plt.ylim(0, 1)
    plt.xlim(ages.min(), ages.max())
    plt.legend(title="Dataset")
    plt.tight_layout()
    plt.show()


def plot_sample_deviations(datasets):
    """
    Generates a ridge density plot to visualize the distribution of sample deviations
    from the population mean for multiple datasets.

    Args:
        datasets (dict): A dictionary of display names to GeoData objects.

    Returns:
        None: Displays the plots directly
    """
    # Ensure inputs are dict for uniform processing
    if not isinstance(datasets, dict):
        raise ValueError("Must pass in valid dictionary")

    # Create a DataFrame to store deviations and dataset IDs
    combined_data = []

    for name, dataset in datasets.items():
        quality_report = dataset.quality_report()
        deviations = quality_report.sample_report["deviation"]
        deviations = deviations.to_frame(name="Deviation")
        deviations["Dataset"] = name
        combined_data.append(deviations)

    # Concatenate all data
    combined_df = pd.concat(combined_data, ignore_index=True)
    # Sort dataset IDs for consistent stacking
    dataset_labels = combined_df["Dataset"].unique()[
        ::-1
    ]  # Reverse order for top-down stacking
    num_datasets = len(dataset_labels)
    # Prepare the figure and axes for subplots
    fig, axes = plt.subplots(
        nrows=num_datasets,
        ncols=1,
        figsize=(10, num_datasets * 1.5),
        sharex=True,
        gridspec_kw={"hspace": 0.5},
    )
    # Loop through each dataset and create density plots
    for idx, (dataset_id, ax) in enumerate(zip(dataset_labels, axes)):
        subset = combined_df[combined_df["Dataset"] == dataset_id]["Deviation"]
        sns.kdeplot(
            subset,
            ax=ax,
            fill=True,
            alpha=0.8,
            linewidth=1,
            color=sns.color_palette("viridis", num_datasets)[idx],
            bw_adjust=0.5,
        )
        ax.set_ylabel(
            dataset_id, fontsize=10, rotation=0, labelpad=40, va="center"
        )
        ax.set_yticks([])
        ax.set_xlim(0, combined_df["Deviation"].max() + 0.01)
        for spine in ["top", "right", "left"]:
            ax.spines[spine].set_visible(False)
    axes[-1].set_xlabel("Deviation", fontsize=12)
    plt.show()


def plot_health_outcome(models, data, health_outcome_col):
    """
    Creates a boxplot and scatterplot to visualize predictions grouped by health outcome.

    Parameters:
    - models: List of models to generate predictions.
    - data: Dataset used for predictions and metadata.
    - health_outcome_col: Column name for the health outcome.

    Returns:
    - None: Displays the plot directly.
    """
    if health_outcome_col is None:
        # Raise an error if health outcome column is not provided
        raise ValueError(
            "health_outcome_col is required for 'health_outcome' plot."
        )

    # Generate predictions for all models
    all_results = get_predictions(models, data)
    # Merge predictions with metadata for visualization
    merged_data = merge_predictions(data.metadata, all_results)

    plt.figure(figsize=(10, 6))
    for model in models:
        name = model.details["name"]
        prediction_col = f"{name}_Predicted"
        # Create a boxplot for predictions grouped by health outcome
        sns.boxplot(
            x=health_outcome_col,
            y=prediction_col,
            data=merged_data,
            hue=health_outcome_col,
            width=0.3,
            legend=False,
        )
        # Overlay scatterplot to show individual data points
        sns.stripplot(
            x=health_outcome_col,
            y=prediction_col,
            data=merged_data,
            hue=health_outcome_col,
            jitter=True,
            dodge=False,
            legend=False,
        )

    plt.ylabel("Predicted Score")
    sns.despine()
    plt.show()


def plot_clock_correlation_matrix(models, data):
    """
    Creates a clustered heatmap of correlations between clock predictions, age, and sex.

    Parameters:
    - models: List of models to generate predictions.
    - data: Dataset used for predictions and metadata.

    Returns:
    - None: Displays the plot directly.
    """
    # Generate predictions for all models
    all_results = get_predictions(models, data)
    merged_data = merge_predictions(data.metadata, all_results)

    # Ensure required columns are present in the dataset
    if "age" not in merged_data.columns or "sex" not in merged_data.columns:
        raise ValueError("The dataset must contain 'age' and 'sex' columns.")

    # Convert 'sex' to numeric codes if it's a categorical variable
    if merged_data["sex"].dtype.name == "category":
        merged_data["sex"] = merged_data["sex"].cat.codes

    # Select columns for correlation analysis
    correlation_columns = [
        col for col in merged_data.columns if col.endswith("_Predicted")
    ] + ["age", "sex"]
    correlation_matrix = merged_data[correlation_columns].corr()

    # Clean up column names for better readability
    cleaned_labels = {
        col: col.replace("_Predicted", "")
        for col in correlation_matrix.columns
    }
    correlation_matrix.rename(
        index=cleaned_labels, columns=cleaned_labels, inplace=True
    )
    # Reorder columns for clustering with 'sex' and 'age' first
    reorder_columns = ["sex", "age"] + [
        col for col in correlation_matrix.columns if col not in ["sex", "age"]
    ]
    reordered_matrix = correlation_matrix.loc[reorder_columns, reorder_columns]

    # Create a clustered heatmap
    sns.clustermap(
        reordered_matrix,
        annot=True,
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        linewidths=0.5,
        fmt=".2f",
        figsize=(12, 12),
        cbar_kws={"label": "Correlation Coefficient"},
        row_cluster=True,
        col_cluster=True,
    )
    plt.show()


def plot_age_prediction(models, data):
    """
    Plots scatterplots for predicted vs actual age.

    Parameters:
    - models: List of models to generate predictions.
    - data: Dataset used for predictions and metadata.

    Returns:
    - None: Displays the plot directly.
    """
    all_results = get_predictions(models, data)
    merged_data = merge_predictions(data.metadata, all_results)

    # Extract prediction columns and ensure 'age' column exists
    prediction_columns = [
        col for col in merged_data.columns if col.endswith("_Predicted")
    ]
    if "age" not in merged_data.columns:
        raise ValueError("The dataset must contain an 'age' column.")

    # Remove rows with missing values in relevant columns
    merged_data = merged_data.dropna(subset=["age"] + prediction_columns)

    plt.figure(figsize=(12, 8))
    for col in prediction_columns:
        # Plot scatterplot for each model's predictions vs actual age
        sns.scatterplot(
            x="age",
            y=col,
            data=merged_data,
            label=col.replace("_Predicted", ""),
            alpha=0.7,
        )
    # Overlay a scatterplot for the actual age as reference
    sns.scatterplot(
        x="age",
        y="age",
        data=merged_data,
        label="Age",
        color="red",
        marker="+",
        s=100,
    )
    plt.xlabel("Chronological Age", fontsize=14)
    plt.ylabel("Predicted Age", fontsize=14)
    plt.legend(title="Models", fontsize=12)
    plt.grid(visible=True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()


def plot_clock_deviation_heatmap(models, data):
    """
    Creates a heatmap of deviations between predicted and chronological age.

    Parameters:
    - models: List of models to generate predictions.
    - data: Dataset used for predictions and metadata.

    Returns:
    - None: Displays the plot directly.
    """
    all_results = get_predictions(models, data)
    merged_data = merge_predictions(data.metadata, all_results)

    # Ensure 'age' column exists in the dataset
    if "age" not in merged_data.columns:
        raise ValueError("The dataset must contain an 'age' column.")

    # Compute deviations for each model's predictions
    prediction_columns = [
        col for col in merged_data.columns if col.endswith("_Predicted")
    ]
    for col in prediction_columns:
        merged_data[f"{col}_Deviation"] = merged_data[col] - merged_data["age"]

    # Extract deviation columns for the heatmap
    deviation_columns = [f"{col}_Deviation" for col in prediction_columns]
    clock_names = [
        col.replace("_Predicted_Deviation", "") for col in deviation_columns
    ]

    # Dynamically calculate figure size based on the number of samples and models
    num_samples = len(merged_data)
    num_models = len(clock_names)
    row_height = 0.2
    col_width = 0.2
    figsize = (
        max(12, col_width * num_models),
        max(12, row_height * num_samples),
    )

    plt.figure(figsize=(figsize))
    sns.heatmap(
        merged_data[deviation_columns],
        cmap="coolwarm",
        center=0,
        annot=True,
        cbar_kws={"label": "Deviation (Epigenetic Age - Chronological Age)"},
        yticklabels=merged_data.index,
        xticklabels=clock_names,
    )
    plt.xlabel("Clocks")
    plt.ylabel("Samples (IDs)")
    plt.tight_layout(pad=2)
    plt.show()


def get_predictions(models, data):
    """
    Generates predictions for a list of models and a dataset.

    Parameters:
    - models: List of models to generate predictions.
    - data: Dataset used for predictions.

    Returns:
    - all_results: List of DataFrames containing predictions from each model.
    """
    all_results = []
    for model in models:
        # Generate predictions for the dataset
        results = model.predict(data)
        name = model.details["name"]
        # Validate that model metadata contains a 'name'
        if "name" not in model.details:
            raise ValueError("Model metadata must contain a 'name'.")
        # Rename prediction column to include the model name
        results.rename(
            columns={"Predicted": f"{name}_Predicted"}, inplace=True
        )
        # Ensure indices are strings for consistent merging
        results.set_index(results.index.astype(str), inplace=True)
        all_results.append(results)
    return all_results


def merge_predictions(metadata, all_results):
    """
    Merges predictions from multiple models with the dataset metadata.

    Parameters:
    - metadata: Metadata DataFrame containing sample information.
    - all_results: List of DataFrames containing model predictions.

    Returns:
    - merged_data: DataFrame with predictions merged into the metadata.
    """
    merged_data = metadata.copy()
    for results in all_results:
        # Merge predictions into metadata
        merged_data = merged_data.join(results, how="left")
    return merged_data
