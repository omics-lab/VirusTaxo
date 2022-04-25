import pickle
from Bio import SeqIO
from collections import defaultdict
import numpy as np


def entropy(x):
    eps = 1e-11
    if len(x) == 1:
        return 0.0

    shift = max(x)
    # print(shift)
    x -= shift
    # print(x)
    summation = sum(np.exp(x))
    # print(summation)
    p = np.exp(x) / summation
    # print(p)
    entropy = sum(p * np.log(p+eps)) / np.log(len(p))
    entropy *= -1
    return entropy


def predict(fasta_file, molecule_type):

    if molecule_type == 'RNA':
        model_path, k = '/path', 17
    elif molecule_type == 'DNA':
        model_path, k = '/path', 21
    else:
        model_path, k = '/path', 20

    # model_path, k = './model/RNA/model_k_17.pkl', 17

    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    print('Id\tLength\tGenus\tEntropy')

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

            if not genus_score_tuple:
                E = 1.0
                return -1, -1, E
            else:
                E = entropy(np.array([x for _, x in genus_score_tuple]))
                prediction = max(genus_score_tuple, key=lambda x: x[1])[0]
                cnt = max(genus_score_tuple, key=lambda x: x[1])[1]
                return prediction, cnt, E

        prediction_1, cnt_1, E_1 = get_genus(record.seq)
        prediction_2, cnt_2, E_2 = get_genus(record.seq.reverse_complement())
        length = len(record.seq)

        if prediction_1 == -1 and prediction_2 == -1:
            print(f'{record.id}\t{length}\tUnclassified\t1.0')
        else:
            if cnt_1 > cnt_2:
                print(f'{record.id}\t{length}\t{prediction_1}\t{E_1}')
            else:
                print(f'{record.id}\t{length}\t{prediction_2}\t{E_2}')


if __name__ == '__main__':
    fasta_file = 'temp.fa'
    predict(fasta_file, 'RNA')
