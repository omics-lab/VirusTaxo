import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from config import params
from tqdm import tqdm
import os
from util import save_object, entropy
import pickle


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', required=True,
                        help='metadata file path')
    parser.add_argument('--seq', required=True,
                        help='Fasta sequence file path')
    parser.add_argument('--model_file', required=True,
                        help='pre saved model path')

    args = parser.parse_args()

    df = pd.read_csv(args.meta)

    accession_to_genus = {row['Accession']: row['Genus'] for _, row in df.iterrows()}

    # db = defaultdict(set)

    with open(args.model_file, 'rb') as f:
        model = pickle.load(f)

    unclassified, classified = 0, 0

    for record in tqdm(SeqIO.parse(args.seq, format='fasta')):
        if record.id not in accession_to_genus:
            continue
        # print(f'predicting for the id: {record.id}')
        seq = str(record.seq)
        count = defaultdict(int)

        # genus = accession_to_genus[record.id]
        for idx in range(len(seq) - params['k'] + 1):
            kmer = seq[idx:idx+params['k']]
            if kmer in model:
                matched_genera = model[kmer]
                if params['discriminative'] and len(matched_genera) == 1:
                    count[list(matched_genera)[0]] += 1
                else:
                    for genus in matched_genera:
                        count[genus] += 1
        genus_score_tuple = list(count.items())

        if not genus_score_tuple:
            E = 1.0
        else:
            E = entropy(np.array([x for _, x in genus_score_tuple]))

        if E > 0.5:
            unclassified += 1
        else:
            prediction = max(genus_score_tuple, key=lambda x: x[1])[0]
            if prediction == accession_to_genus[record.id]:
                classified += 1

    total = len(accession_to_genus)
    # print(unclassified)
    # print(classified)
    # print("N\tUnclassified\tRemaining\tCorrectly_Classified\tAccuracy")
    n = total - unclassified
    print(f'{total}\t{unclassified}\t{n}\t{classified}\t{classified / n * 100}')

