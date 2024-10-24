from biolearn.model import model_definitions
from biolearn.util import get_data_file
import yaml
import csv
import os

def generate_models_csv(models, output_file="generated/model_table.csv"):
    header = ["Name", "Year", "Species", "Tissue", "Predicts", "Source", "Coefficients"]

    # Generate the CSV data
    with open(output_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header
        writer.writerow(header)

        # Add each model to the CSV
        for name, details in models.items():
            if details["source"] == "unknown":
                source_link = "unknown"
            else:
                source_link = f"`paper <{details['source']}>`_"
            coefficients_link = f"`coefficients file <https://github.com/bio-learn/biolearn/blob/master/biolearn/data/{details['model']['file']}>`_"

            row = [
                name,
                details["year"],
                details["species"],
                details["tissue"],
                details.get('output', 'unknown'),
                source_link,
                coefficients_link
            ]
            writer.writerow(row)

    print(f"CSV generated at: {output_file}")

def ensure_folder_exists(path):
    # Get the directory part of the path (if it is a file)
    dir_path = os.path.dirname(path)

    if os.path.isdir(path):
        dir_path = path

    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def generate_data_csv_from_yaml(yaml_file, output_file="generated/data_table.csv"):
    header = ["ID", "Title", "Format", "Samples", "Age Present", "Sex Present"]

    with open(yaml_file, 'r') as stream:
        yaml_data = yaml.safe_load(stream)

    with open(output_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)

        for item in yaml_data["items"]:
            age_present = "Yes" if "age" in item.get("parser", {}).get("metadata", {}) else "No"
            sex_present = "Yes" if "sex" in item.get("parser", {}).get("metadata", {}) else "No"
            geo_link = f"`{item['id']} <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={item['id']}>`_"
            title = item['title'][:60] + '...' if len(item['title']) > 60 else item['title']
            row = [
                geo_link,
                title,
                item["format"],
                item["samples"],
                age_present,
                sex_present
            ]
            writer.writerow(row)

    print(f"CSV generated at: {output_file}")

def generate_model_usage_csv(models, output_file="generated/model_usage.csv"):
    # Define the CSV header
    header = ["Name", "Commercial Usage", "Non-Commercial Usage"]

    # Ensure the folder exists
    ensure_folder_exists(output_file)

    # Generate the CSV data
    with open(output_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header
        writer.writerow(header)

        # Add each model's usage data to the CSV
        for name, details in models.items():
            usage = details.get("usage", {})
            source = details.get("source", "unknown")

            # Handle commercial usage
            commercial_usage = usage.get("commercial", "unknown")
            if commercial_usage == "unknown" and source != "unknown":
                commercial_usage = f"Contact `paper <{source}>`_ author"
            
            # Handle non-commercial usage
            non_commercial_usage = usage.get("non-commercial", "unknown")
            if non_commercial_usage == "unknown":
                non_commercial_usage = "Free to use"

            row = [
                name,
                commercial_usage,
                non_commercial_usage
            ]
            writer.writerow(row)

    print(f"Model usage CSV generated at: {output_file}")



if __name__ == "__main__":
    ensure_folder_exists("generated/")
    generate_models_csv(model_definitions)
    generate_model_usage_csv(model_definitions)
    generate_data_csv_from_yaml(get_data_file("library.yaml"))
