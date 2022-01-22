import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from config import params
from tqdm import tqdm
import os
from util import save_object


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', required=True,
                        help='metadata file path')
    parser.add_argument('--seq', required=True,
                        help='Fasta sequence file path')
    parser.add_argument('--model_file', required=True,
                        help='path for saving model')

    args = parser.parse_args()

    df = pd.read_csv(args.meta)

    accession_to_genus = {row['Accession']: row['Genus'] for _, row in df.iterrows()}

    db = defaultdict(set)

    for record in tqdm(SeqIO.parse(args.seq, format='fasta')):
        if record.id not in accession_to_genus:
            continue
        seq = str(record.seq)
        genus = accession_to_genus[record.id]
        for idx in range(len(seq) - params['k'] + 1):
            db[seq[idx:idx+params['k']]].add(genus)

    dirname = os.path.dirname(args.model_file)

    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    save_object(db, args.model_file)

    # print(db)







