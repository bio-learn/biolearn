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

class MethylationAnalyzer:
    def __init__(self, datasets, dataset_ids):
        """
        Initializes the MethylationAnalyzer with multiple datasets.

        Parameters:
        - datasets (list of tuples): Each tuple contains two DataFrames, (data_dnam, data_meta).
            * data_dnam: DNAm data with CpG sites as rows and samples as columns.
            * data_meta: Metadata with sample IDs as index, containing 'age' and 'sex' columns.
        - dataset_ids (list of str): List of dataset IDs corresponding to the datasets.
        """
        self.datasets = datasets
        self.dataset_ids = dataset_ids  # Store dataset IDs for proper labeling

    def identify_stable_cpg_sites(self, threshold=0.01):
        """
        Identifies stable CpG sites in each dataset based on variance.

        Parameters:
        - threshold: Variance threshold to determine stability (default: 0.01)

        Returns:
        - stable_sites_summary: List of dictionaries with stable CpG sites per dataset.
        """
        stable_sites_summary = []
        for (data_dnam, _), dataset_id in zip(self.datasets, self.dataset_ids):
            variances = data_dnam.var(axis=1)
            stable_sites = variances[variances < threshold].index.tolist()
            stable_sites_summary.append({'Dataset': dataset_id, 'Stable CpG Sites': stable_sites})
        return stable_sites_summary

    def compute_age_group_stats(self, cpg_site, age_bins=5):
        """
        Computes statistics (mean, median, std, etc.) for methylation levels by age group and sex.

        Parameters:
        - cpg_site: CpG site ID to compute statistics for.
        - age_bins: Number of age bins for grouping.

        Returns:
        - stats_df: DataFrame containing computed statistics.
        - combined_data: DataFrame with added 'age_group' column.
        """
        combined_data = []
        sex_mapping = {1: 'M', 2: 'F', '1': 'M', '2': 'F', 'M': 'M', 'F': 'F', 'Male': 'M', 'Female': 'F'}

        # Combine data from all datasets
        for (data_dnam, data_meta), dataset_id in zip(self.datasets, self.dataset_ids):
            if cpg_site not in data_dnam.index:
                print(f"Warning: CpG site '{cpg_site}' not found in dataset {dataset_id}. Skipping this dataset.")
                continue
            methylation_data = data_dnam.loc[cpg_site]
            valid_samples = data_meta.index.intersection(methylation_data.index)
            plot_data = data_meta.loc[valid_samples, ['age', 'sex']].copy()
            plot_data['methylation'] = methylation_data[valid_samples]
            plot_data['sex'] = plot_data['sex'].replace(sex_mapping)
            plot_data['dataset'] = dataset_id
            combined_data.append(plot_data)

        combined_data = pd.concat(combined_data).dropna(subset=['age', 'methylation'])

        # Define bins for age groups with whole numbers
        min_age = int(np.floor(combined_data['age'].min()))
        max_age = int(np.ceil(combined_data['age'].max()))
        bins = np.linspace(min_age, max_age, age_bins + 1).astype(int)
        combined_data['age_group'] = pd.cut(
            combined_data['age'], 
            bins=bins, 
            include_lowest=True, 
            labels=[f"{bins[i]}-{bins[i+1]-1}" for i in range(len(bins)-1)]
        )

        stats_summary = []

        # Group by age_group to compute statistics for each age range
        for age_group in combined_data['age_group'].unique():
            group_data = combined_data[combined_data['age_group'] == age_group]
            
            # Separate data by sex within the age group
            male_data = group_data[group_data['sex'] == 'M']['methylation']
            female_data = group_data[group_data['sex'] == 'F']['methylation']
            
            # Calculate statistics
            stats_row = {
                'age_group': age_group,
                'male_mean': male_data.mean(),
                'female_mean': female_data.mean(),
                'male_median': male_data.median(),
                'female_median': female_data.median(),
                'male_std': male_data.std(),
                'female_std': female_data.std(),
                'male_count': male_data.count(),
                'female_count': female_data.count()
            }

            # Append statistics for the current age group
            stats_summary.append(stats_row)

        # Convert the summary list to a DataFrame for easy viewing
        stats_df = pd.DataFrame(stats_summary)
        return stats_df, combined_data

    def plot_methylation_by_age_sex(self, combined_data, cpg_site):
        """
        Visualizes methylation levels at a specific CpG site by age group and sex using combined data.

        Parameters:
        - combined_data: DataFrame with combined data from all datasets, including 'age_group' column.
        - cpg_site: CpG site ID to visualize.
        """
        plt.figure(figsize=(14, 8))
        sns.violinplot(x='age_group', y='methylation', hue='sex', data=combined_data, split=True, palette="Set2")
        plt.xlabel('Age Group')
        plt.ylabel(f'Methylation Level at {cpg_site}')
        plt.legend(title="Sex")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()

    def compute_methylation_stats(self, cpg_site, degree=2):
        """
        Computes regression statistics (linear and polynomial) for DNA methylation vs. age.

        Parameters:
        - cpg_site: CpG site ID to compute statistics for.
        - degree: Degree of the polynomial for the polynomial fit (default: 2).

        Returns:
        - combined_data: DataFrame with combined methylation and age data.
        - statistics: Dictionary containing linear and polynomial regression results.
        """
        combined_data = []

        for (data_dnam, data_meta), dataset_id in zip(self.datasets, self.dataset_ids):
            if cpg_site not in data_dnam.index:
                print(f"Warning: CpG site '{cpg_site}' not found in dataset {dataset_id}. Skipping this dataset.")
                continue

            # Extract methylation data for the CpG site
            methylation_data = data_dnam.loc[cpg_site]
            valid_samples = data_meta.index.intersection(methylation_data.index)
            plot_data = data_meta.loc[valid_samples, ['age']].copy()
            plot_data['methylation'] = methylation_data[valid_samples]
            plot_data['dataset'] = dataset_id  # Use dataset ID instead of generic name
            combined_data.append(plot_data)

        # Combine data across datasets
        combined_data = pd.concat(combined_data).dropna(subset=['age', 'methylation'])

        # Prepare data for regression
        ages = combined_data['age'].values.reshape(-1, 1)
        methylation_levels = combined_data['methylation'].values

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
            'linear_pred': linear_pred,
            'r_squared_linear': r_squared_linear,
            'corr': corr,
            'p_value': p_value,
            'poly_pred': poly_pred,
            'r_squared_poly': r_squared_poly,
            'degree': degree
        }

        # Print statistics summary
        print(f"Computed statistics for CpG site: {cpg_site}")
        print(f"Linear Regression R^2: {r_squared_linear:.4f}")
        print(f"Pearson Correlation (r): {corr:.4f}")
        print(f"Pearson Correlation p-value: {p_value:.4e}")
        print(f"Polynomial Regression (Degree {degree}) R^2: {r_squared_poly:.4f}")

        return combined_data, statistics

    def plot_methylation_vs_age(self, combined_data, statistics, cpg_site):
        """
        Plots the DNA methylation level of a specific CpG site against age,
        with both a linear and polynomial regression line.

        Parameters:
        - combined_data: DataFrame with combined methylation and age data.
        - statistics: Dictionary containing regression statistics.
        - cpg_site: CpG site ID to visualize.
        """
        ages = combined_data['age'].values
        methylation_levels = combined_data['methylation'].values

        plt.figure(figsize=(12, 6))

        # Scatter plot of methylation levels, color-coded by dataset
        sns.scatterplot(
            x=ages,
            y=methylation_levels,
            hue=combined_data['dataset'],
            palette="tab10",
            s=10,
            alpha=0.6,
            legend='full'
        )

        # Linear fit line
        sns.lineplot(x=ages, y=statistics['linear_pred'], color='red', linestyle="--")

        # Polynomial fit line
        sns.lineplot(x=ages, y=statistics['poly_pred'], color='green')

        plt.xlabel("Age")
        plt.ylabel("DNA Methylation")
        plt.ylim(0, 1)
        plt.xlim(ages.min(), ages.max())
        plt.legend(title="Dataset")
        plt.tight_layout()
        plt.show()


class ModelAnalyzer:
    def __init__(self):
        self.gallery = ModelGallery()

    def plot_health_outcome(self, models, data, health_outcome_col):
        # Raise an error if the health_outcome_col is not provided
        if health_outcome_col is None:
            raise ValueError("health_outcome_col is required for 'health_outcome' plot.")

        # Generate predictions for the models and merge with the dataset metadata
        all_results = self._get_predictions(models, data)
        merged_data = self._merge_predictions(data.metadata, all_results)

        # Create a boxplot and scatterplot to visualize predictions grouped by health outcome
        plt.figure(figsize=(10, 6))
        for model in models:
            name = model.metadata["name"]
            prediction_col = f"{name}_Predicted"
            sns.boxplot(
                x=health_outcome_col,
                y=prediction_col,
                data=merged_data,
                hue=health_outcome_col,
                width=0.3,
                legend=False
            )
            sns.stripplot(
                x=health_outcome_col,
                y=prediction_col,
                data=merged_data,
                hue=health_outcome_col,
                jitter=True,
                dodge=False,
                legend=False
            )
        
        # Add labels and finalize the plot
        plt.ylabel("Predicted Score")
        sns.despine()
        plt.show()

    def plot_clock_correlation_matrix(self, models, data):
        all_results = self._get_predictions(models, data)
        merged_data = self._merge_predictions(data.metadata, all_results)

        # Ensure that necessary columns ('age' and 'sex') are present in the dataset
        if "age" not in merged_data.columns or "sex" not in merged_data.columns:
            raise ValueError("The dataset must contain 'age' and 'sex' columns.")

        # Convert 'sex' column to numeric if it's categorical
        if merged_data["sex"].dtype.name == "category":
            merged_data["sex"] = merged_data["sex"].cat.codes

        # Identify columns for correlation analysis and calculate the correlation matrix
        correlation_columns = [col for col in merged_data.columns if col.endswith("_Predicted")] + ["age", "sex"]
        correlation_matrix = merged_data[correlation_columns].corr()
        
        # Clean column labels for better readability
        cleaned_labels = {col: col.replace("_Predicted", "") for col in correlation_matrix.columns}
        correlation_matrix.rename(index=cleaned_labels, columns=cleaned_labels, inplace=True)
        reorder_columns = ['sex', 'age'] + [col for col in correlation_matrix.columns if col not in ['sex', 'age']]
        reordered_matrix = correlation_matrix.loc[reorder_columns, reorder_columns]

        # Plot the clustered heatmap
        sns.clustermap(
            reordered_matrix,
            annot=True, 
            cmap="coolwarm",  
            vmin=-1, vmax=1,  
            linewidths=0.5,  
            fmt=".2f",  
            figsize=(12, 12), 
            cbar_kws={"label": "Correlation Coefficient"},  
            row_cluster=True, 
            col_cluster=True,  
        )
        plt.show()

    def plot_age_prediction(self, models, data):
        all_results = self._get_predictions(models, data)
        merged_data = self._merge_predictions(data.metadata, all_results)

        # Identify prediction columns and ensure the 'age' column is present
        prediction_columns = [col for col in merged_data.columns if col.endswith("_Predicted")]
        if "age" not in merged_data.columns:
            raise ValueError("The dataset must contain an 'age' column.")

        # Drop rows with missing values in relevant columns
        merged_data = merged_data.dropna(subset=["age"] + prediction_columns)

        # Plot scatterplots for predicted vs actual age
        plt.figure(figsize=(12, 8))
        for col in prediction_columns:
            sns.scatterplot(
                x="age",
                y=col,
                data=merged_data,
                label=col.replace("_Predicted", ""),
                alpha=0.7,
            )
        sns.scatterplot(
            x="age",
            y="age",
            data=merged_data,
            label="Age",
            color="red",
            marker="+",
            s=100
        )
        plt.xlabel("Chronological Age", fontsize=14)
        plt.ylabel("Predicted Age", fontsize=14)
        plt.legend(title="Models", fontsize=12)
        plt.grid(visible=True, linestyle="--", alpha=0.6)
        plt.tight_layout()
        plt.show()

    def plot_clock_deviation_heatmap(self, models, data):
        all_results = self._get_predictions(models, data)
        merged_data = self._merge_predictions(data.metadata, all_results)

        # Ensure the 'age' column is present
        if "age" not in merged_data.columns:
            raise ValueError("The dataset must contain an 'age' column.")

        # Calculate deviations for each model's predictions
        prediction_columns = [col for col in merged_data.columns if col.endswith("_Predicted")]
        for col in prediction_columns:
            merged_data[f"{col}_Deviation"] = merged_data[col] - merged_data["age"]

        # Extract deviation columns and format for heatmap
        deviation_columns = [f"{col}_Deviation" for col in prediction_columns]
        clock_names = [col.replace("_Predicted_Deviation", "") for col in deviation_columns]

        # Plot heatmap
        plt.figure(figsize=(10, len(merged_data) * 0.2))
        sns.heatmap(
            merged_data[deviation_columns],
            cmap="coolwarm",
            center=0,
            annot=False,
            cbar_kws={"label": "Deviation (Epigenetic Age - Chronological Age)"},
            yticklabels=merged_data.index,
            xticklabels=clock_names,  
        )
        plt.xlabel("Epigenetic Clocks")
        plt.ylabel("Samples (IDs)")
        plt.tight_layout(pad=2)
        plt.show()

    def _get_predictions(self, models, data):
        all_results = []
        for model in models:
            results = model.predict(data)
            name = model.metadata["name"]
            if "name" not in model.metadata:
                raise ValueError("Model metadata must contain a 'name'.")
            # Rename prediction column to include model name and set index as string
            results.rename(columns={"Predicted": f"{name}_Predicted"}, inplace=True)
            results.set_index(results.index.astype(str), inplace=True)
            all_results.append(results)
        return all_results

    def _merge_predictions(self, metadata, all_results):
        merged_data = metadata.copy()
        for results in all_results:
            # Join predictions with metadata
            merged_data = merged_data.join(results, how="left")
        return merged_data

