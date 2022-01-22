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
    df = df[df['Genus'].notna()]

    mp = defaultdict(list)

    for _, row in df.iterrows():
        # print(dict(row))
        mp[row['Genus']].append(dict(row))

    train, test = [], []

    for key, val in mp.items():
        seq_num = len(val)
        if seq_num < 3:
            continue
        temp = val
        np.random.shuffle(temp)
        n = max(1, int(seq_num * 0.1))

        test.extend(temp[0:n])
        train.extend(temp[n:])

    # print(len(train), len(test))
    train_df = pd.DataFrame(train)
    test_df = pd.DataFrame(test)
    print(train_df.shape)
    print(test_df.shape)


    # train_df.to_csv(os.path.join(args.out_dir, 'train.csv'), header=True, index=False)
    # test_df.to_csv(os.path.join(args.out_dir, 'test.csv'), header=True, index=False)





    # print(df.groupby(by=["Genus"])['Genus'].count())