import pickle
from Bio import SeqIO
from collections import defaultdict
import numpy as np
from util import entropy
import argparse


def predict(fasta_file, model_path):
    # print(f'Loading the model...')

    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    print('Id\tLength\tGenus\tEntropy\tEnrichment_Score')
    # get k-mer length
    for key in model.keys():
        k = len(key)
        break

    for record in SeqIO.parse(fasta_file, "fasta"):
        def get_genus(seq):
            count = defaultdict(int)
            for idx in range(len(seq) - k + 1):
                kmer = seq[idx:idx + k]
                if kmer in model:
                    matched_genera = model[kmer]
                    for genus in matched_genera:
                        count[genus] += 1
            genus_score_tuple = list(count.items())
            # print(f'PrintTouple\t{genus_score_tuple}') 
            total_kmers = len(seq) - k + 1 
            # print(total_kmers)

            if not genus_score_tuple:
                E = 1.0 # no genera are found
                es = 0 # no genera are found
                return -1, -1, E, es
            else:
                E = entropy(np.array([x for _, x in genus_score_tuple]))
                prediction = max(genus_score_tuple, key=lambda x: x[1])[0]
                cnt = max(genus_score_tuple, key=lambda x: x[1])[1]
                es = cnt / total_kmers
                return prediction, cnt, E, es

        prediction_1, cnt_1, E_1, es1 = get_genus(record.seq)
        prediction_2, cnt_2, E_2, es2 = get_genus(record.seq.reverse_complement())
        length = len(record.seq)

        if prediction_1 == -1 and prediction_2 == -1:
            print(f'{record.id}\t{length}\tNoHit\t1.0\t0')
        else:
            if cnt_1 > cnt_2:
                print(f'{record.id}\t{length}\t{prediction_1}\t{E_1:.3f}\t{es1:.3f}')
            else:
                print(f'{record.id}\t{length}\t{prediction_2}\t{E_2:.3f}\t{es2:.3f}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--model_path', required=True,
                        help='Absolute or relative path of pre-built model')

    parser.add_argument('--seq', required=True,
                        help='Absolute or relative path of fasta sequence file')

    args = parser.parse_args()

    predict(args.seq, args.model_path)