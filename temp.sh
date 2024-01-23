# test

python3 build.py \
   --meta ./Dataset/RNA_meta.csv \
   --seq ./Dataset/RNA_seq.fasta \
   --k 17 \
   --saving_path ./RNA.pkl

python3 predict.py \
   --model_path ./RNA.pkl \
   --seq ./input.fasta

# DNA
python3 build.py \
   --meta ./Dataset/DNA_9384_meta.csv \
   --seq ./Dataset/sequences_20240122_3701960.fasta \
   --k 21 \
   --saving_path ./DNA_9384_k21.pkl

# RNA
python3 build.py \
   --meta ./Dataset/DNA_9384_meta.csv \
   --seq ./Dataset/sequences_20240122_3701960.fasta \
   --k 21 \
   --saving_path ./DNA_9384_k21.pkl 
   
# DNA_RNA
python3 build.py \
   --meta ./Dataset/DNA_9384_meta.csv \
   --seq ./Dataset/sequences_20240122_3701960.fasta \
   --k 21 \
   --saving_path ./DNA_9384_k21.pkl 