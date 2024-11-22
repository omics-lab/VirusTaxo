import pickle
import numpy as np
from collections import defaultdict


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def get_rank(seq, model, k):
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


if __name__ == '__main__':
    E = entropy(np.array([1, 0, 0, 6, 1]))
    print(E)