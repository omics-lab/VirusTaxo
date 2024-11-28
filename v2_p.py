import pickle
from Bio import SeqIO
import pandas as pd
import csv
from util import get_rank
from tqdm import tqdm
import os
import argparse


def predict(fasta_file, database_path, output_csv, entropy, enrichment, enrichment_spp):
    # Define rank names and initialize list to store output file paths
    ranks = ["Family", "Genus", "Species"]
    temp_dir = "./temp/"
    os.makedirs(temp_dir, exist_ok=True)
    output_files = []

    for rank in ranks:
        # Set the default enrichment based on rank
        if rank == "Species":
            enrichment_value = enrichment_spp  # Use species-specific enrichment if provided
        else:
            enrichment_value = enrichment  # Use the general enrichment value for Genus and Family

        # Print the entropy and enrichment values for the current rank
        print(f"Processing rank: {rank}")
        print(f"Using entropy threshold: {entropy}")
        print(f"Using enrichment threshold: {enrichment_value}")

        model_file = os.path.join(database_path, f"{rank}_database.pkl")
        output_file = os.path.join(temp_dir, f"{rank}_predictions.csv")
        output_files.append(output_file)

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
                    if E_1 >= entropy or es1 <= enrichment_value:
                        prediction_1 = "Unclassified"
                    if E_2 >= entropy or es2 <= enrichment_value:
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

    # Step 1: Load input files
    # Load metadata
    metadata_file = os.path.join(database_path, "metadata.csv")
    metadata = pd.read_csv(metadata_file)
    metadata = metadata[["Family", "Genus", "Species"]]

    # Step 2: Create sets for each row of metadata
    metadata_sets = metadata.apply(lambda row: set(row), axis=1).tolist()

    # Step 3: Create sets for each row of merged_df
    merged_df['RowSet'] = merged_df[["Family", "Genus", "Species"]].apply(lambda row: set(row), axis=1)

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

    print(f"Number of rows with valid taxonomy: {yes_count}")
    print(f"Number of rows with invalid taxonomy: {no_count}")

    # Drop the temporary 'RowSet' column and save the output
    merged_df.drop(columns=['RowSet'], inplace=True)
    merged_df.to_csv(output_csv, index=False)

    print(f"Updated file saved to: {output_csv}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--database_path', required=True,
                        help='Absolute or relative path containing three model files (Family_database.pkl, Genus_database.pkl and Species_database.pkl)')

    parser.add_argument('--seq', required=True,
                        help='Absolute or relative path of fasta sequence file')

    parser.add_argument('--output_csv', required=False, default="VirusTaxo_merged_predictions.csv",
                        help='Path to save the merged output CSV file (default: VirusTaxo_predictions.csv)')

    parser.add_argument('--entropy', required=False, default=0.5, type=float,
                        help='Entropy threshold (default: 0.5)')

    parser.add_argument('--enrichment', required=False, default=0.05, type=float,
                        help='Enrichment score threshold for Genus and Family (default: 0.05)')

    parser.add_argument('--enrichment_spp', required=False, default=0.8, type=float,
                        help='Enrichment score threshold for Species (default: 0.8)')

    args = parser.parse_args()

    predict(args.seq, args.database_path, args.output_csv, args.entropy, args.enrichment, args.enrichment_spp)
