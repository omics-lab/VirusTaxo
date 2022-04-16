import pickle
import numpy as np


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


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