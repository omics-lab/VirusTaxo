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