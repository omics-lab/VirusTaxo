import random
import pandas as pd
from Bio import SeqIO

def get_accessions_from_metadata(metadata_file):
    # Read the metadata CSV file and extract the Accession column
    metadata = pd.read_csv(metadata_file)
    return set(metadata["Accession"].values)

def split_fasta(input_fasta, train_fasta, test_fasta, metadata_file, train_percentage=0.8):
    # Get the set of accession IDs from the metadata file
    accession_set = get_accessions_from_metadata(metadata_file)
    
    # Read all sequences from the input FASTA file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Filter sequences that are present in the accession set
    filtered_sequences = [seq for seq in sequences if seq.id in accession_set]
    
    # Shuffle the filtered sequences randomly
    random.shuffle(filtered_sequences)
    
    # Calculate the split index based on the train_percentage
    split_index = int(len(filtered_sequences) * train_percentage)
    
    # Split the sequences into training and testing sets
    train_sequences = filtered_sequences[:split_index]
    test_sequences = filtered_sequences[split_index:]
    
    # Write the training sequences to the train_fasta file
    with open(train_fasta, "w") as train_file:
        SeqIO.write(train_sequences, train_file, "fasta")
    
    # Write the testing sequences to the test_fasta file
    with open(test_fasta, "w") as test_file:
        SeqIO.write(test_sequences, test_file, "fasta")
    
    print(f"FASTA file has been split into {train_fasta} (80%) and {test_fasta} (20%).")

# Input and output file paths
input_fasta = "./temp/sequences.fasta"
train_fasta = "./temp/train.fasta"
test_fasta = "./temp/test.fasta"
metadata_file = "./temp/metadata.csv"  # Path to the metadata CSV file

# Split the FASTA file
split_fasta(input_fasta, train_fasta, test_fasta, metadata_file)
