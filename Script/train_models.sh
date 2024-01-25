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
   --meta ./Dataset/DNA_RNA_18451_meta.csv \
   --seq ./Dataset/sequences_20240122_3701960.fasta \
   --k 20 \
   --saving_path ./DNA_RNA_18451_k20.pkl 

# RNA
python3 build.py \
   --meta ./Dataset/RNA_meta_4k.csv \
   --seq ./Dataset/RNA_seq.fasta \
   --k 17 \
   --saving_path ./RNA_4k.pkl






# 
python3 predict.py \
    --model_path ./DNA_9384_k21.pkl \
    --seq ./Dataset/input.fasta

Id      Length  Genus   Entropy
NC_004205.1 (RNA virus)    280     Lazarusvirus    0.0 (DNA virus)
NC_038276.1     280     Unclassified    1.0

# 
python3 predict.py \
    --model_path ./RNA_9067_k17.pkl \
    --seq ./Dataset/input.fasta >out.txt 

Id      Length  Genus   Entropy
NC_004205.1     280     Seadornavirus   0.0
NC_038276.1     280     Ledantevirus    0.0

# 
python3 predict.py \
    --model_path ./DNA_RNA_18451_k20.pkl \
    --seq ./Dataset/input.fasta
#accurate 
Id      Length  Genus   Entropy
NC_004205.1     280     Seadornavirus   0.0
NC_038276.1     280     Ledantevirus    0.0

#
python3 predict.py \
    --model_path ./DNA_RNA_18451_k20.pkl \
    --seq ./Dataset/asm_islam.et.al_6648_covid.fasta >asm_islam.et.al_6648_covid.fasta.out.txt  

#
python3 predict.py \
    --model_path ./DNA_9384_k21.pkl \
    --seq ./Dataset/asm_islam.et.al_6648_covid.fasta >asm_islam.et.al_6648_covid.fasta.DNAmodel.out.txt      

/projects/epigenomics3/temp/rislam/VirusTaxo/