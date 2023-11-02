from biolearn.model import clock_definitions
import csv
import os

def generate_clocks_csv(clocks, output_file="generated/clock_table.csv"):
    header = ["Name", "Year", "Species", "Tissue", "Source", "Coefficients"]

    # Generate the CSV data
    with open(output_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header
        writer.writerow(header)

        # Add each clock to the CSV
        for name, details in clocks.items():
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


if __name__ == "__main__":
    ensure_folder_exists("generated")
    generate_clocks_csv(clock_definitions)