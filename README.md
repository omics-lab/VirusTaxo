### Installation

```
git clone https://github.com/omics-lab/VirusTaxo_v1
cd VirusTaxo_v1
python3 -m venv environment
source ./environment/bin/activate
pip install -r requirements.txt
```

### Validation

```
mkdir -p model/DNA
bash cross_validate.sh
```

### Final RNA or DNA
```
mkdir -p ./model/Final_RNA
python3 select_data.py --meta ./Dataset/RNA/RNA_meta.csv --out_path ./Dataset/RNA/meta_fi.tsv
python3 train.py --meta ./Dataset/RNA/meta_fi.tsv --seq ./Dataset/RNA/RNA_seq.fasta --model_dir ./model/Final_RNA
```
### Metagenomics

```
python3 benchmark.py --megahit_dir ./megahit --model_file ./model/Final_RNA/model_k_17.pkl
```


### Train Test Split

```
python3 train_test_spiltter.py --meta ./Dataset/RNA/RNA_meta.csv \
    --out_dir ./Dataset/RNA
```
### Build the DB
```
 python3 train.py --meta ./Dataset/RNA/train.csv --seq ./Dataset/RNA/RNA_seq.fasta --model_dir ./model/RNA
```

### Accuracy Find

```
python3 acc_find.py --meta ./Dataset/RNA/test.csv \
    --seq ./Dataset/RNA/RNA_seq.fasta \
    --model_file ./model/RNA/model.pkl
```