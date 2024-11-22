import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from tqdm import tqdm
from util import save_object


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', required=True,
                        help='File path of metadata')
    parser.add_argument('--seq', required=True,
                        help='Fasta sequence file path')
    parser.add_argument('--k', required=True, type=int,
                        help='Length of k-mer')
    parser.add_argument('--saving_path', required=True,
                        help='Absolute or relative path for saving model')

    args = parser.parse_args()

    k = args.k

    df = pd.read_csv(args.meta)

    accession_to_genus = {row['Id']: row['Genus'] for _, row in df.iterrows()}

    db = defaultdict(set)

    for record in tqdm(SeqIO.parse(args.seq, format='fasta')):
        if record.id not in accession_to_genus:
            continue
        seq = str(record.seq)
        genus = accession_to_genus[record.id]
        for idx in range(len(seq) - k + 1):
            db[seq[idx : idx + k]].add(genus)

    # print('Length before discard:', len(db))
    discard_candidate = set()

    for key, val in db.items():
        if len(val) > 1:
            discard_candidate.add(key)

    for key in discard_candidate:
        db.pop(key)

    save_object(db, args.saving_path)
