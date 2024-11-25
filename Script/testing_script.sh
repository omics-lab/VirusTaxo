python3 predict.py \
   --model_path ../temp/vt_db_jan21_2024/RNA*.pkl \
   --seq ./Dataset/asm_islam.et.al_6648_covid_random.fasta >../temp/testout

python3 predict.py \
	--model_path ../../temp/vt_db_jan21_2024/DNA_RNA*pkl \
	--seq ../Dataset/asm_head.fasta 

less ../temp/testout | head 
less ../temp/testout | awk '$4 <=.5 && $5 >= .8{print $0}' | head 


python v2_b.py \
	 --meta ./Dataset/Accession_Species_Genus_Family_12612_meta_1k.csv \
	 --seq ../temp/vt2_database/sequences_20240122_3701960.fasta \
	 --k 20 \
	 --saving_dir ../temp/vt2_database

python3 v2_p.py \
	--database_path ../temp/vt2_database/ \
	--seq ./Dataset/asm_head.fasta 

# check acc python at 5-fold 80% vs 20%

# check with different kmers


#!/bin/bash

## k-mer acc

# cd /path/VirusTaxo/

python split_fasta.py

# train.fasta, test.fasta 

for k in {15..17}
do
  echo "Current k-mer: $k";
  
  # Run the first Python script
  python v2_b.py \
    --meta ./temp/metadata.csv \
    --seq ./temp/seq1k.fasta \
    --k "$k" \
    --saving_dir ./temp/k_mer_loop/

  echo "Predicting with k-mer = $k";

  # Run the second Python script
  python3 v2_p.py \
    --database_path ./temp/k_mer_loop/ \
    --seq ./temp/seq100.fasta \
    --output_csv ./temp/k_mer_loop/VirusTaxo_merged_predictions_kmer_"$k".csv
done
