import pickle
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import csv
from util import entropy, get_rank
from tqdm import tqdm
import argparse


def predict(fasta_file, model_path, output_csv="VirusTaxo_taxonomy_output.csv", custom_E=0.5, custom_es=0.8):
    # Load the model
    print("Loading the model...")
    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    # Count the number of sequences in the FASTA file
    num_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
    print(f"Total number of query sequences = {num_sequences}")

    # Open the CSV file for writing
    with open(output_csv, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write the header
        csv_writer.writerow(['Id', 'Length', 'Genus', 'Entropy', 'Enrichment_Score'])

        # Get k-mer length
        for key in model.keys():
            k = len(key)
            break
        print(f"k-mer size = {k}")

        # Process each record in the FASTA file with a progress bar
        with tqdm(total=num_sequences, desc="Processing sequences") as pbar:
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

    print(f"Predictions saved to {output_csv}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--model_path', required=True,
                        help='Absolute or relative path of pre-built model')

    parser.add_argument('--seq', required=True,
                        help='Absolute or relative path of fasta sequence file')

    parser.add_argument('--output_csv', required=False, default="VirusTaxo_taxonomy_output.csv",
                        help='Path to save the output CSV file (default: VirusTaxo_taxonomy_output.csv)')

    parser.add_argument('--entropy', required=False, type=float, default=0.5,
                        help='Entropy threshold for classification (default: 0.5)')

    parser.add_argument('--enrichment_score', required=False, type=float, default=0.8,
                        help='Enrichment score threshold for classification (default: 0.8)')

    args = parser.parse_args()

    predict(args.seq, args.model_path, args.output_csv, args.entropy, args.enrichment_score)
