#!/usr/bin/env bash

for i in {1..5}; do
  echo "Iteration: ${i}"
  echo "spliting the dataset"
  python3 80_20_splitter.py --meta ./Dataset/DNA/DNA_meta.csv --out_dir ./Dataset/DNA
  echo "training..."
  python3 train.py --meta ./Dataset/DNA/train.csv --seq ./Dataset/DNA/DNA_seq.fasta --model_dir ./model/DNA
  echo "testing.."
  python3 acc_find.py --meta ./Dataset/DNA/test.csv --seq ./Dataset/DNA/DNA_seq.fasta \
    --model_file ./model/DNA/model_k_17.pkl >> DNA_cross.tsv
done