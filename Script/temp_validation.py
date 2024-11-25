import pandas as pd
import os

# Load the datasets
merged_df_file = os.path.join("../Dataset/", "test_merged_df.csv")
merged_df = pd.read_csv(merged_df_file)

metadata_file = os.path.join("../Dataset/", "test_metadata.csv")
metadata = pd.read_csv(metadata_file)
metadata = metadata[["Species", "Genus", "Family"]]  # Ensure only relevant columns

# Convert metadata into a set of tuples for fast lookup
metadata_set = set(metadata.apply(tuple, axis=1))

def is_valid_taxonomy(row, metadata_set):
    species = row["Species"]
    genus = row["Genus"]
    family = row["Family"]
    
    # Check if all three ranks exist in metadata
    if (species, genus, family) in metadata_set:
        return f"('{species}', '{genus}', '{family}')"

    # Check if one rank is 'Unclassified' and the other two form a valid pair
    if species == "Unclassified" and (genus, family) in metadata_set:
        return f"('{species}', '{genus}', '{family}')"
    if genus == "Unclassified" and (species, family) in metadata_set:
        return f"('{species}', '{genus}', '{family}')"
    if family == "Unclassified" and (species, genus) in metadata_set:
        return f"('{species}', '{genus}', '{family}')"

    # Check if two ranks are 'Unclassified' and the third is valid
    if species == "Unclassified" and genus == "Unclassified" and family in metadata_set:
        return f"('{species}', '{genus}', '{family}')"
    if genus == "Unclassified" and family == "Unclassified" and species in metadata_set:
        return f"('{species}', '{genus}', '{family}')"
    if species == "Unclassified" and family == "Unclassified" and genus in metadata_set:
        return f"('{species}', '{genus}', '{family}')"

    # Check if all three ranks are 'Unclassified'
    if species == "Unclassified" and genus == "Unclassified" and family == "Unclassified":
        return f"('{species}', '{genus}', '{family}')"

    # If no valid taxonomy, return 'Wrong classification'
    return "('Wrong classification')"

# Apply the function to the merged_df and store the result in the 'Valid_Taxonomy' column
merged_df["Valid_Taxonomy"] = merged_df.apply(lambda row: is_valid_taxonomy(row, metadata_set), axis=1)

# Save the result to a new CSV file
merged_df.to_csv("output_with_valid_taxonomy.csv", index=False)
