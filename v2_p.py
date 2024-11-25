import pickle
from Bio import SeqIO
import pandas as pd
import csv
from util import get_rank
from tqdm import tqdm
import os
import argparse


def predict(fasta_file, database_path, output_csv="VirusTaxo_merged_predictions.csv", custom_E=0.5, custom_es=0.8):
    # Define rank names and initialize list to store output file paths
    ranks = ["Species", "Genus", "Family"]
    temp_dir = "./temp/"
    os.makedirs(temp_dir, exist_ok=True)
    output_files = []

    # Load metadata
    metadata_file = os.path.join(database_path, "metadata.csv")
    metadata = pd.read_csv(metadata_file)
    metadata = metadata[["Species", "Genus", "Family"]]  # Ensure only relevant columns
    metadata_set = set(metadata.apply(tuple, axis=1))  # Create a set of valid (Species, Genus, Family) combinations

    for rank in ranks:
        model_file = os.path.join(database_path, f"{rank}_database.pkl")
        output_file = os.path.join(temp_dir, f"{rank}_predictions.csv")
        output_files.append(output_file)

        print(f"Processing taxonomic rank: {rank}")
        print(f"Loading database from: {model_file}")

        # Load the model
        with open(model_file, 'rb') as f:
            model = pickle.load(f)

        # Count the number of sequences in the FASTA file
        num_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        print(f"Total number of query sequences loaded for {rank} prediction = {num_sequences}")

        # Open the CSV file for writing
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            # Write the header
            csv_writer.writerow(['Accession', 'Query_Seq_Length', rank, f'{rank}_Entropy', f'{rank}_Enrichment'])

            # Get k-mer length
            for key in model.keys():
                k = len(key)
                break
            print(f"k-mer size for {rank} prediction = {k}")

            # Process each record in the FASTA file with a progress bar
            with tqdm(total=num_sequences, desc=f"Processing query sequences for {rank}") as pbar:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    # Get predictions for both the sequence and its reverse complement
                    prediction_1, cnt_1, E_1, es1 = get_rank(record.seq, model, k)
                    prediction_2, cnt_2, E_2, es2 = get_rank(record.seq.reverse_complement(), model, k)
                    length = len(record.seq)

                    # Check if the conditions for "Unclassified" are met for both predictions using custom thresholds
                    if E_1 >= custom_E or es1 <= custom_es:
                        prediction_1 = "Unclassified"
                    if E_2 >= custom_E or es2 <= custom_es:
                        prediction_2 = "Unclassified"

                    # Write the results based on the best prediction
                    if prediction_1 == -1 and prediction_2 == -1:
                        csv_writer.writerow([record.id, length, 'NoHit', 1.0, 0])
                    else:
                        if cnt_1 > cnt_2:
                            csv_writer.writerow([record.id, length, prediction_1, f'{E_1:.3f}', f'{es1:.3f}'])
                        else:
                            csv_writer.writerow([record.id, length, prediction_2, f'{E_2:.3f}', f'{es2:.3f}'])

                    pbar.update(1)  # Update the progress bar for each record

        print(f"Predictions saved to {output_file}")

    # Merge the output CSV files
    print("Merging predictions from all ranks...")
    merged_df = None
    for i, rank_file in enumerate(output_files):
        df = pd.read_csv(rank_file)
        if i == 0:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on=["Accession", "Query_Seq_Length"], how="outer")

    # Save merged predictions to CSV before filtering
    merged_csv_file = os.path.join(temp_dir, "merged_predictions_before_filtering.csv")
    merged_df.to_csv(merged_csv_file, index=False)
    print(f"Merged predictions saved to {merged_csv_file} before filtering.")

    # Validate and filter merged predictions based on taxonomic hierarchy
    print("Validating and filtering merged predictions based on valid taxonomic hierarchy...")

    # Function to check valid taxonomy for non-unclassified ranks
    # def is_valid_taxonomy(row, metadata_set):
    #     species = row["Species"]
    #     genus = row["Genus"]
    #     family = row["Family"]

    #     # Check if all three ranks exist in metadata
    #     if (species, genus, family) in metadata_set:
    #         return (species, genus, family)

    #     # Check for valid pairs when one rank is 'Unclassified'
    #     if species == "Unclassified" and (genus, family) in metadata_set:
    #         return (species, genus, family)
    #     if genus == "Unclassified" and (species, family) in metadata_set:
    #         return (species, genus, family)
    #     if family == "Unclassified" and (species, genus) in metadata_set:
    #         return (species, genus, family)

    #     # Additional checks for when two or three ranks are 'Unclassified'
    #     if species == "Unclassified" and genus == "Unclassified" and family in metadata_set:
    #         return (species, genus, family)
    #     if genus == "Unclassified" and family == "Unclassified" and species in metadata_set:
    #         return (species, genus, family)
    #     if species == "Unclassified" and family == "Unclassified" and genus in metadata_set:
    #         return (species, genus, family)

    #     if species == "Unclassified" and genus == "Unclassified" and family == "Unclassified":
    #         return (species, genus, family)

    #     return None

    def is_valid_taxonomy(row, metadata_set):
        species = row["Species"]
        genus = row["Genus"]
        family = row["Family"]

        print(f"Checking row: {species}, {genus}, {family}")

        # Check if all three ranks exist in metadata
        if (species, genus, family) in metadata_set:
            print(f"Found valid taxonomy: {species}, {genus}, {family}")
            return (species, genus, family)

        # Check for valid pairs when one rank is 'Unclassified'
        if species == "Unclassified" and (genus, family) in metadata_set:
            print(f"Found valid pair for species: {species}, {genus}, {family}")
            return (species, genus, family)
        if genus == "Unclassified" and (species, family) in metadata_set:
            print(f"Found valid pair for genus: {species}, {genus}, {family}")
            return (species, genus, family)
        if family == "Unclassified" and (species, genus) in metadata_set:
            print(f"Found valid pair for family: {species}, {genus}, {family}")
            return (species, genus, family)

        # Check if two ranks are 'Unclassified'
        if species == "Unclassified" and genus == "Unclassified" and family in metadata_set:
            print(f"Found valid pair for species/genus: {species}, {genus}, {family}")
            return (species, genus, family)
        if genus == "Unclassified" and family == "Unclassified" and species in metadata_set:
            print(f"Found valid pair for genus/family: {species}, {genus}, {family}")
            return (species, genus, family)
        if species == "Unclassified" and family == "Unclassified" and genus in metadata_set:
            print(f"Found valid pair for species/family: {species}, {genus}, {family}")
            return (species, genus, family)

        # Check if all three ranks are 'Unclassified'
        if species == "Unclassified" and genus == "Unclassified" and family == "Unclassified":
            print(f"All ranks are Unclassified: {species}, {genus}, {family}")
            return (species, genus, family)

        return None

    merged_df["Valid_Taxonomy"] = merged_df.apply(lambda row: is_valid_taxonomy(row, metadata_set), axis=1)
    merged_df.to_csv("output_with_valid_taxonomy.csv", index=False)

    # # Apply the validation function
    # merged_df["Valid_Taxonomy"] = merged_df.apply(lambda row: is_valid_taxonomy(row, metadata_set), axis=1)

    # initial_row_count = len(merged_df)  # Count initial number of rows
    # filtered_df = merged_df[merged_df["Valid_Taxonomy"]].copy()  # Only keep rows with valid taxonomy
    # filtered_row_count = len(filtered_df)  # Count number of rows after filtering
    # rows_filtered_out = initial_row_count - filtered_row_count

    # print(f"Ambiguous predictions are filtered out: {rows_filtered_out}")
    # print(f"Successful predictions of query sequences: {filtered_row_count}")
    # filtered_df.drop(columns=["Valid_Taxonomy"], inplace=True)  # Drop the temporary column

    # # Save the filtered merged CSV file
    # filtered_df.to_csv(output_csv, index=False)
    # print(f"Filtered and merged predictions saved to {output_csv}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--database_path', required=True,
                        help='Absolute or relative path containing model files (e.g., Species_database.pkl, Genus_database.pkl, Family_database.pkl)')

    parser.add_argument('--seq', required=True,
                        help='Absolute or relative path of fasta sequence file')

    parser.add_argument('--output_csv', required=False, default="VirusTaxo_merged_predictions.csv",
                        help='Path to save the merged output CSV file (default: VirusTaxo_merged_predictions.csv)')

    parser.add_argument('--custom_E', required=False, default=0.5, type=float,
                        help='Custom entropy threshold (default: 0.5)')

    parser.add_argument('--custom_es', required=False, default=0.8, type=float,
                        help='Custom enrichment score threshold (default: 0.8)')

    args = parser.parse_args()

    predict(args.seq, args.database_path, args.output_csv, args.custom_E, args.custom_es)
