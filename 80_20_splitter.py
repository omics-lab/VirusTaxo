import argparse
import pandas as pd
from collections import defaultdict
import numpy as np
import datetime as dt
import os


np.random.seed(dt.datetime.now().microsecond)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', required=True,
                        help='Raw metadata file path')
    parser.add_argument('--out_dir', required=True,
                        help='output dir path for train.csv & test.csv')

    args = parser.parse_args()

    df = pd.read_csv(args.meta)
    # df = df.head(100)
    # We Skip if genus information is missing
    df = df[df['Genus'].notna() & df['Species'].notna()]

    print(df.shape)

    mp = defaultdict(int)

    for _, row in df.iterrows():
        mp[row['Genus']] += 1

    discard = [key for key, val in mp.items() if val < 3]

    print(len(discard))

    for val in discard:
        df = df[df.Genus != val]

    print(df.shape)

    list_of_dict = df.to_dict('records')

    np.random.shuffle(list_of_dict)

    train_n = int(len(list_of_dict) * 0.8)

    train = pd.DataFrame(list_of_dict[0:train_n])
    test = pd.DataFrame(list_of_dict[train_n:])

    train.to_csv(os.path.join(args.out_dir, 'train.csv'), header=True, index=False)
    test.to_csv(os.path.join(args.out_dir, 'test.csv'), header=True, index=False)
