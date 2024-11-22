import os
import pickle
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import argparse
import mmap

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def entropy(x):
    eps = 1e-11
    if len(x) == 1:
        return 0.0

    shift = max(x)
    x -= shift
    summation = sum(np.exp(x))
    p = np.exp(x) / summation
    entropy = sum(p * np.log(p + eps)) / np.log(len(p))
    entropy *= -1
    return entropy

def get_genus(seq, model, k, return_counts=False):
    count = defaultdict(int)
    for idx in range(len(seq) - k + 1):
        kmer = seq[idx:idx + k]
        if kmer in model:
            matched_genera = model[kmer]
            for genus in matched_genera:
                count[genus] += 1

    genus_score_tuple = list(count.items())
    total_kmers = len(seq) - k + 1

    if not genus_score_tuple:
        return -1, -1, 1.0, 1.0 if not return_counts else count
    else:
        E = entropy(np.array([x for _, x in genus_score_tuple]))
        prediction = max(genus_score_tuple, key=lambda x: x[1])[0]
        cnt = max(genus_score_tuple, key=lambda x: x[1])[1]
        es = cnt / total_kmers
        return (prediction, cnt, E, es) if not return_counts else count

def split_model(model_path, output_dir, chunk_size):
    """
    Splits a large model file into smaller chunks.

    Parameters:
    - model_path (str): Path to the input .pkl file.
    - output_dir (str): Directory to save the chunks.
    - chunk_size (int): Number of entries per chunk.
    """
    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    os.makedirs(output_dir, exist_ok=True)
    items = list(model.items())
    total_chunks = (len(items) + chunk_size - 1) // chunk_size

    for i in range(total_chunks):
        chunk = dict(items[i * chunk_size:(i + 1) * chunk_size])
        chunk_file = os.path.join(output_dir, f'model_chunk_{i + 1}.pkl')
        with open(chunk_file, 'wb') as f:
            pickle.dump(chunk, f, pickle.HIGHEST_PROTOCOL)

    print(f"Model split into {total_chunks} chunks and saved to '{output_dir}'.")

def predict_with_chunks(fasta_file, chunk_dir):
    """
    Runs predictions by iterating over all model chunks.

    Parameters:
    - fasta_file (str): Path to the input fasta sequence file.
    - chunk_dir (str): Directory containing the model chunks.
    """
    chunk_files = sorted([os.path.join(chunk_dir, f) for f in os.listdir(chunk_dir) if f.endswith('.pkl')])
    if not chunk_files:
        raise ValueError("No model chunks found in the specified directory.")

    print('Id\tLength\tGenus\tEntropy\tEnrichment_Score')

    # Process each fasta record
    for record in SeqIO.parse(fasta_file, "fasta"):
        length = len(record.seq)
        aggregated_counts = defaultdict(int)

        # Process each chunk using memory-mapped model
        for chunk_file in chunk_files:
            # Use mmap to open the chunk in a memory-mapped way
            with open(chunk_file, 'rb') as f:
                mmapped_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                # Read and load the chunk
                model = pickle.loads(mmapped_file[:])

            # Get k-mer length from the model
            k = len(next(iter(model.keys())))

            # Update aggregated counts for forward and reverse sequences
            forward_counts = get_genus(record.seq, model, k, return_counts=True)
            reverse_counts = get_genus(record.seq.reverse_complement(), model, k, return_counts=True)

            for genus, count in forward_counts.items():
                aggregated_counts[genus] += count
            for genus, count in reverse_counts.items():
                aggregated_counts[genus] += count

        # Final prediction based on aggregated counts
        genus_score_tuple = list(aggregated_counts.items())
        if not genus_score_tuple:
            print(f'{record.id}\t{length}\tNoHit\t1.0\t0')
        else:
            E = entropy(np.array([x for _, x in genus_score_tuple]))
            prediction = max(genus_score_tuple, key=lambda x: x[1])[0]
            cnt = max(genus_score_tuple, key=lambda x: x[1])[1]
            total_kmers = len(record.seq) - k + 1
            es = cnt / total_kmers
            print(f'{record.id}\t{length}\t{prediction}\t{E:.3f}\t{es:.3f}')


def main():
    parser = argparse.ArgumentParser(description="Split model and predict using chunks with memory-mapping.")

    # Add arguments
    parser.add_argument('--model_path', required=True, help='Path to the input .pkl model file.')
    parser.add_argument('--seq', required=True, help='Path to the fasta sequence file.')
    parser.add_argument('--chunk_dir', required=True, help='Directory to store model chunks or containing existing chunks.')
    parser.add_argument('--chunk_size', type=int, default=100000, help='Number of entries per chunk (only used if splitting).')

    # Parse arguments
    args = parser.parse_args()

    # Split the model into chunks if necessary, and then predict
    if not os.path.exists(args.chunk_dir) or not any(f.endswith('.pkl') for f in os.listdir(args.chunk_dir)):
        print("Splitting the model into chunks...")
        split_model(args.model_path, args.chunk_dir, args.chunk_size)

    print("Running predictions using the model chunks...")
    predict_with_chunks(args.seq, args.chunk_dir)


if __name__ == '__main__':
    main()
