import os
import pandas as pd

# Step 1: Load input files
metadata_file = os.path.join("../Dataset/", "test_metadata.csv")
metadata = pd.read_csv(metadata_file)
metadata = metadata[["Species", "Genus", "Family"]]

merged_df_file = os.path.join("../Dataset/", "test_merged_df.csv")
merged_df = pd.read_csv(merged_df_file)

# Step 2: Create sets for each row of metadata
metadata_sets = metadata.apply(lambda row: set(row), axis=1).tolist()

# Step 3: Create sets for each row of merged_df
merged_df['RowSet'] = merged_df[["Species", "Genus", "Family"]].apply(lambda row: set(row), axis=1)

# Step 4: Define a function to determine validity strictly per row
def is_valid(row_set):
    for metadata_set in metadata_sets:
        unclassified_count = list(row_set).count("Unclassified")
        # Check if the row matches a metadata row or partially matches based on conditions
        if unclassified_count == 3:  # All are "Unclassified"
            return "Yes"
        elif unclassified_count == 2:  # Two are "Unclassified"
            return "Yes"
        elif unclassified_count == 1:  # One is "Unclassified"
            classified_elements = row_set - {"Unclassified"}
            if classified_elements.issubset(metadata_set):
                return "Yes"
        elif row_set == metadata_set:  # Exact match
            return "Yes"
    return "No"

# Apply validity function to merged_df
merged_df['Valid'] = merged_df['RowSet'].apply(is_valid)

yes_count = (merged_df['Valid'] == "Yes").sum()
no_count = (merged_df['Valid'] == "No").sum()

print(f"Number of rows with 'Yes': {yes_count}")
print(f"Number of rows with 'No': {no_count}")

# Drop the temporary 'RowSet' column and save the output
merged_df.drop(columns=['RowSet'], inplace=True)
output_file = os.path.join("../Dataset/", "test_merged_df_with_valid.csv")
merged_df.to_csv(output_file, index=False)

print(f"Updated file saved to: {output_file}")
