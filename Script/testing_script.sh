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


# build final model

python v2_b.py \
  --meta ./temp/database/metadata.csv \
  --seq ./temp/database/sequences.fasta \
  --k 16 \
  --saving_dir ./temp/database/

# check acc python at 5-fold 80% vs 20%

# check with different kmers


#!/bin/bash

## k-mer acc

# cd /path/VirusTaxo/

python split_fasta.py

# train.fasta, test.fasta 

# test
for k in {10..14}
do
  echo "Current k-mer: $k";
  
  # Run the first Python script
  python v2_b.py \
    --meta ../temp/metadata.csv \
    --seq ../temp/seq1k.fasta \
    --k "$k" \
    --saving_dir ../temp/;

	echo "Predicting with k-mer = $k";

  # Run the second Python script
  python3 v2_p.py \
    --database_path ../temp/ \
    --seq ../temp/seq100.fasta \
    --output_csv VirusTaxo_predictions_"$k".csv;
done

  python3 v2_p.py \
    --database_path ../temp/ \
    --seq ../temp/seq100.fasta \
    --output_csv VirusTaxo_predictions.csv;

# in server 
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/VirusTaxo/
#scp -r /home/rashedul/project/VirusTaxo/ rislam@gphost03.bcgsc.ca:/projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/

source activate py310
python -m venv environment
source ./environment/bin/activate
pip install -r requirements.txt

# 80% 20%
python split_fasta.py

# train.fasta, test.fasta 

for k in {10..14}
do
  echo "Current k-mer: $k";
  
  # Run the first Python script
  python v2_b.py \
    --meta ./temp/metadata.csv \
    --seq ./temp/train.fasta \
    --k "$k" \
    --saving_dir ./temp/k_mer_loop/;

	echo "Predicting with k-mer = $k";

  # Run the second Python script
  python3 v2_p.py \
    --database_path ./temp/k_mer_loop/ \
    --seq ./temp/test.fasta \
    --output_csv ./temp/k_mer_loop/VirusTaxo_predictions_"$k".csv;
done >vt_kmers_.5.05._k10-14.log


# 3-fold cross validation

for k in {1..3}
do
	echo $k;
	python split_fasta.py;

  # Run the first Python script
  python v2_b.py \
    --meta ./temp/metadata.csv \
    --seq ./temp/train.fasta \
    --k 16 \
    --saving_dir ./temp/cv/;

  # Run the second Python script
  python3 v2_p.py \
    --database_path ./temp/cv/ \
    --seq ./temp/test.fasta \
    --output_csv ./temp/cv/VirusTaxo_predictions_fold_"$k".csv;
done >vt_3cv.log

# acc testing
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/VirusTaxo/temp/k_mer_loop/

# invalid 
cat VirusTaxo_predictions_*.csv | grep -v Yes | wc 

# count classified
for file in VirusTaxo_predictions*.csv; do
    echo $file
    for f in 3 6 9; do
        echo "rank: $f:"
        less $file | awk -v col=$f -F, '{print $col}' | grep -v Uncl | wc -l
    done
done | paste - - - - - - - 


# acc test0
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/VirusTaxo/temp/k_mer_loop/

# fam acc
for f in VirusTaxo_predictions_*.csv;
do 
	echo $f;
	less $f | awk -F, '{print $1"_"$3}' | sort >temp1
	a=$(less temp1 | grep -v Uncl | wc -l);
	echo $a;
	less ../metadata.csv | tr -d '"' | awk -F, '{print $2"_"$5}' | sort >temp2
	b=$(grep -f temp1 temp2 | wc -l);
	echo $b;
	acc=$(awk "BEGIN {print $b / $a}");
	echo "accuracy:" $acc;
	rm temp1 temp2 
done | paste - - - - 

# gen acc
for f in VirusTaxo_predictions_*.csv;
do 
	echo $f;
	less $f | awk -F, '{print $1"_"$6}' | sort >temp1
	a=$(less temp1 | grep -v Uncl | wc -l);
	echo $a;
	less ../metadata.csv | tr -d '"' | awk -F, '{print $2"_"$4}' | sort >temp2
	b=$(grep -f temp1 temp2 | wc -l);
	echo $b;
	acc=$(awk "BEGIN {print $b / $a}");
	echo "accuracy:" $acc;
	rm temp1 temp2 
done | paste - - - - 

# spp acc
for f in VirusTaxo_predictions_*.csv;
do 
	echo $f;
	less $f | awk -F, '{print $1"_"$9}' | sort >temp1
	a=$(less temp1 | grep -v Uncl | wc -l);
	echo $a;
	less ../metadata.csv | tr -d '"' | awk -F, '{print $2"_"$3}' | sort >temp2
	b=$(grep -f temp1 temp2 | wc -l);
	echo $b;
	acc=$(awk "BEGIN {print $b / $a}");
	echo "accuracy:" $acc;
	rm temp1 temp2 
done | paste - - - - 
