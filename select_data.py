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
    parser.add_argument('--out_path', required=True,
                        help='output path for final meta')

    args = parser.parse_args()

    df = pd.read_csv(args.meta)
    # df = df.head(100)
    # We Skip if genus information is missing
    df = df[df['Genus'].notna() & df['Species'].notna()]

    # print(df.shape)

    mp = defaultdict(int)

    for _, row in df.iterrows():
        mp[row['Genus']] += 1

    discard = [key for key, val in mp.items() if val < 3]

    # print(len(discard))

    for val in discard:
        df = df[df.Genus != val]

    # print(df.shape)
    df.to_csv(args.out_path, header=True, index=False)

