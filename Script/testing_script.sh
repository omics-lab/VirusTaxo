#### build final model
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/VirusTaxo/

python build.py \
  --meta ./temp/database/metadata.csv \
  --seq ./temp/database/sequences.fasta \
  --k 16 \
  --saving_dir ./temp/database/

#### final model validation
# asm test was successful: all genus and family and 80% species detection with default
# random seqs were not classified at default
# no other taxonomy was predicted in both covid and random fasta
# in the metavirome data, 43% sequences had predicted genus with ES=0 but with default ES it was 1%. so ES cutoff removed significant number of seqs. 
# in genecode 1k lncRNA, 45% sequences had predicted genus with ES=0 but with default ES it was 0.3%. so ES cutoff removed significant number of seqs. 

for f in ./temp/test/*meta*fasta;
do echo $f;
python3 v2_p.py \
   --database_path ./temp/database/ \
   --seq $f \
   --output_csv $f.csv
done

#
for f in ./temp/test/gencode.1k.fasta;
do echo $f;
python3 v2_p.py \
   --database_path ./temp/database/ \
   --seq $f \
   --enrichment 0 \
   --enrichment_spp 0 \
   --output_csv $f.csv
done

#### 3-fold cross validation

for k in {1..2}
do
	echo $k;
	python split_fasta.py;

  # Run the first Python script
  python v2_b_10.py \
    --meta ./temp/metadata.csv \
    --seq ./temp/train.fasta \
    --k 16 \
    --saving_dir ./temp/cv_b10/;

  # Run the second Python script
  python3 v2_p.py \
    --database_path ./temp/cv_b10/ \
    --seq ./temp/test.fasta \
    --output_csv ./temp/cv_b10/VirusTaxo_predictions_fold_10_"$k".csv;
done >vt_3cv_b10.log

#### build models
python v2_b_10.py \
  --meta ./temp/database10/metadata.csv \
  --seq ./temp/database10/sequences.fasta \
  --k 16 \
  --saving_dir ./temp/database10/

#### check acc at 5-fold 80% vs 20%


#!/bin/bash

## k-mer acc

# cd /path/VirusTaxo/

python split_fasta.py

# train.fasta, test.fasta 

# test on local pc

for k in {16..16}
do
  echo "Current k-mer: $k";
  
  # Run the first Python script
  python build.py \
    --meta ../temp/metadata.csv \
    --seq ../temp/seq1k.fasta \
    --k "$k" \
    --saving_dir ../temp/;

	echo "Predicting with k-mer = $k";

  # Run the second Python script
  python3 predict.py \
    --database_path ../temp/ \
    --seq ../temp/seq100.fasta \
    --output_csv ./temp/VirusTaxo_predictions_"$k".csv;
done

python3 v2_p.py \
    --database_path ../temp/ \
    --seq ../temp/seq100.fasta \
    --output_csv VirusTaxo_predictions_16.csv;

#### in server 
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


#### 3-fold cross validation

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
done >vt_3cv_2.log

#### acc testing
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/VirusTaxo/temp/k_mer_loop/

# invalid 
cat VirusTaxo_predictions_*.csv | grep -v Yes | wc 

# count classified
for file in *.csv; do
    echo $file
    for f in 3 6 9; do
        echo "rank: $f:"
        less $file | awk -v col=$f -F, '{print $col}' | grep -v Uncl | wc -l
    done
done | paste - - - - - - - 


# acc test0
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/VirusTaxo/temp/k_mer_loop/
metadata=/projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/VirusTaxo/temp/metadata.csv

# fam acc
for f in *.csv;
do 
	echo $f;
	less $f | awk -F, '{print $1"_"$3}' | sort >temp1
	a=$(less temp1 | grep -v Uncl | wc -l);
	echo $a;
	less $metadata | tr -d '"' | awk -F, '{print $2"_"$5}' | sort >temp2
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
