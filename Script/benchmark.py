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

    parser.add_argument('--megahit_dir', required=True,
                        help='megahit directory file path')

    parser.add_argument('--model_file', required=True,
                        help='pre saved model path')

    args = parser.parse_args()

    with open(args.model_file, 'rb') as f:
        model = pickle.load(f)

    for root, dirs, files in os.walk(args.megahit_dir):
        for dir in dirs:

            fasta_file = os.path.join(root, dir, dir + '.contigs.fa')
            # print(fasta_file)
            for record in SeqIO.parse(fasta_file, "fasta"):
                def get_genus(seq):
                    count = defaultdict(int)
                    for idx in range(len(seq) - params['k'] + 1):
                        kmer = seq[idx:idx + params['k']]
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
                        return -1, -1, E
                    else:
                        E = entropy(np.array([x for _, x in genus_score_tuple]))
                        prediction = max(genus_score_tuple, key=lambda x: x[1])[0]
                        cnt = max(genus_score_tuple, key=lambda x: x[1])[1]
                        return prediction, cnt, E

                seq = str(record.seq)
                rc = record.seq.reverse_complement()
                # print(f'{record.id}')
                prediction_1, cnt_1, E_1 = get_genus(seq)
                prediction_2, cnt_2, E_2 = get_genus(rc)
                if prediction_1 == -1 and prediction_2 == -1:
                    print(f'{dir}\t{record.id}\t{len(seq)}\tUnclassified\t0\t1.0')
                else:
                    if cnt_1 > cnt_2:
                        print(f'{dir}\t{record.id}\t{len(seq)}\t{prediction_1}\t{cnt_1}\t{E_1}')
                    else:
                        print(f'{dir}\t{record.id}\t{len(seq)}\t{prediction_2}\t{cnt_2}\t{E_2}')
