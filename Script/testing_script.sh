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

# task: write python script
I have two inputs:

test_metadata.csv

"Accession,Species,Genus,Family
x1,a,b,c
x2,d,e,f"

test_merged_df.csv

"Entropy,Species,Genus,Family
2,a,b,c
3,Unclassified,b,c
4,a,Unclassified,c
5,a,b,Unclassified
6,Unclassified,Unclassified,Unclassified
7,d,b,c
8,d,Unclassified,c
9,a,b,e"

step-1: load input files from thsoe paths

metadata_file = os.path.join("../Dataset/", "test_metadata.csv")
metadata = pd.read_csv(metadata_file)
metadata = metadata[["Species", "Genus", "Family"]]

merged_df_file = os.path.join("../Dataset/", "test_merged_df.csv")
merged_df = pd.read_csv(merged_df_file)

Step-2: take the three columns Species,Genus,Family from test_metadata.csv and create set for each row of metadata - such (a,b,c) and (d,e,f)

Step-3: similarly take the three columns Species,Genus,Family from test_merged_df.csv and make set for each row for merged_df such (a,b,c) and (Unclassified,b,c)

Step-4: write a new column on test_merged_df.csv named "Valid" and write "Yes" or "No":

write "Yes", if
- merged_df's (Species,Genus,Family) set intersect with metadata's (Species,Genus,Family) set  
- merged_df contains "Unclassified" in one of the (Species,Genus,Family) set and rest of the two elements intersect with metadata's (Species,Genus,Family) set 
- merged_df contains "Unclassified" in two of the (Species,Genus,Family) set  
- merged_df contains three "Unclassified" in (Species,Genus,Family) set  
- otherwise, print "No"

