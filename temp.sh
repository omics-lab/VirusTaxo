#!/bin/bash

# DNA
python3 build.py \
   --meta ./Dataset/DNA_9384_meta.csv \
   --seq ./Dataset/sequences_20240122_3701960.fasta \
   --k 21 \
   --saving_path ./DNA_9384_k21.pkl

# RNA
python3 build.py \
   --meta ./Dataset/RNA_9067_meta.csv \
   --seq ./Dataset/sequences_20240122_3701960.fasta \
   --k 17 \
   --saving_path ./RNA_9067_k17.pkl 
   
# DNA_RNA
python3 build.py \
   --meta ./Dataset/DNA_9384_meta.csv \
   --seq ./Dataset/sequences_20240122_3701960.fasta \
   --k 20 \
   --saving_path ./DNA_RNA_18451_k20.pkl 