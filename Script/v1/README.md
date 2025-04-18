## VirusTaxo: Taxonomic classification of viruses from metagenomic contigs

VirusTaxo has an average accuracy of 93% at genus level across DNA and RNA viruses.

### 1. Running VirusTaxo 
#### Requirements 
- python >= 3.8
- Linux

#### Installation
 - Cloning the repository
```
git clone https://github.com/omics-lab/VirusTaxo
```
 - Creation of Python Virtual Environment
```
cd VirusTaxo
python3 -m venv environment
source ./environment/bin/activate
```
 - Installation of Python Packages
```
pip install -r requirements.txt
```

### 2. Predict virus taxonomy from fasta file using prebuilt database

- Download prebuilt databse (`database.v2_2024`) of VirusTaxo

```
gdown "https://drive.google.com/uc?id=1gz0n5oHomWjpT0HXsrqh8hTLqmqiqgJs"

# Extract db files
tar -xvzf vt_db_jan21_2024.tar.gz
```

- Database files and usage:

| db file | Molecule | Usage |
|----------|----------|----------|
| DNA_RNA_18451_k20.pkl  | DNA & RNA  | Recommended for samples containing both DNA & RNA viruses   |
| DNA_9384_k21.pkl  | DNA  | Recommended for samples containing DNA viruses only |
| RNA_9067_k17.pkl  | RNA  | Recommended for samples containing RNA viruses only |


- Assemble the metagenomic contigs from your metavirome or metagenomic library. 
   - Perform *de novo* assembly using [MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884): 

   ```
   megahit -1 file_R1.fq -2 file_R2.fq --min-contig-len 500 -o contig.fasta
   ```

- Usage of `predict.py`:

```
python3 predict.py -h

usage: predict.py [-h] --model_path MODEL_PATH --seq SEQ [--output_csv OUTPUT_CSV] [--entropy ENTROPY] [--enrichment_score ENRICHMENT_SCORE]

options:
  -h, --help            show this help message and exit
  --model_path MODEL_PATH
                        Absolute or relative path of pre-built model
  --seq SEQ             Absolute or relative path of fasta sequence file
  --output_csv OUTPUT_CSV
                        Path to save the output CSV file (default: VirusTaxo_taxonomy_output.csv)
  --entropy ENTROPY     Entropy threshold for classification (default: 0.5)
  --enrichment_score ENRICHMENT_SCORE
                        Enrichment score threshold for classification (default: 0.8)
```

- Run with the sample [contig.fasta](./Dataset/contig.fasta) file

```
python3 predict.py \
   --model_path /path/DNA_RNA_18451_k20.pkl \ # database file
   --seq ./Dataset/contig.fasta # query fasta file 
```

- Example output for 4 query sequences

```
Id              Length  Genus           Entropy Enrichment_Score
QuerySeq-1      219     Unclassified    1.000   0.000
QuerySeq-2      720     Betacoronavirus 0.000   0.973
QuerySeq-3      1540    Unknown         0.285   0.820
QuerySeq-4      1330    Lentivirus      0.000   0.987
```

### 3. Interpretation of output
- In the taxonomic rank column 
   - `NoHit`: no k-mer overlap between the query and database.
   - `Unknown`: genus name is not assigned in the [ICTV classification](https://ictv.global/). 
   - `Unclassified`: `Entropy` (default >= 0.5) or `Enrichment_Score` (default <= 0.8) is outside of cutoff.

- Lower `Entropy` (such as ≤=0.5) provides the higher level of prediction certainty. You can decrease `Entropy` cutoff for better prediction. 

- Higher `Enrichment_Score` (such as >= 0.8) provides the higher level of prediction certainty. You can increase `Enrichment_Score` cutoff for better prediction. `Enrichment_Score` is the total number of k-mers mapped to the genera divided by total number of k-mers in the query sequence.

### 4. Build your custom database

- Preparing a metadata file in `csv` format. The metadata file must contain columns named `Id`  and `Genus`. Example of metadata file is [here](./Dataset/RNA_meta.csv):

- The sequnce IDs must match with the metadata IDs. Example of input fasta file is [here](./Dataset/RNA_seq.fasta).

 - Building database:

```
python3 build.py \
   --meta ./Dataset/RNA_meta.csv \ # provide your metadata file
   --seq ./Dataset/RNA_seq.fasta \ # provide your fasta file
   --k 17 \
   --saving_path /path/RNA.pkl
```

 - Parameters 
  
   - `meta`: Absolute or relative path of metadata file.
   - `seq`: Absolute or relative path of fasta sequence file.
   - `k` : It denotes the length of k-mer.
   - `saving_path`: Path to save a pickle file (A database File).


### 5. Method limitation and interpretation

- VirusTaxo is trained on known virus sequences and designed to predict taxonomy of virus sequences. 

- Since VirusTaxo uses k-mer enrichment, non-virus sequences could be classified as virus due to random k-mer match. To use VirusTaxo, please make sure to remove non-virus sequences.    

- If your sample contains non-virus sequences, it is highly recommeded to filter out non-viral sequences using tools like [blast](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/find-data/sequence) or [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) `dvf.py -i contig.fasta -o ./`. 

- To avoid contaimination with host sequences, please filter out those by mapping the reads to host reference genomes before using VirusTaxo. 

### 6. Database versions

| Script | Version  | Date     | Sequences | Download |
|----------|----------|----------|----------|----------|
| [v1](Script/v1) | database.v2_2024  | Jan21_2024  | DNA=9384 &  RNA=9067  | [here](https://drive.google.com/file/d/1gz0n5oHomWjpT0HXsrqh8hTLqmqiqgJs/view?usp=sharing)  |
| [v1](Script/v1) | database.v1_2022  | Apr27_2022  | DNA=4421 &  RNA=2529  | [here](https://drive.google.com/file/d/1j9rcFi6AMjA7tSqSizAQO7GpZw-brauZ/view?usp=sharing)  |
| [Used in manuscript](https://github.com/omics-lab/VirusTaxo_Hierarchical) | database.v1_2022  | Apr27_2022  | DNA=4421 &  RNA=2529  | [here](https://drive.google.com/file/d/1j9rcFi6AMjA7tSqSizAQO7GpZw-brauZ/view?usp=sharing) |

### 7. Contact
Rashedul Islam, PhD (rashedul.gen@gmail.com)

### 8. Citation

Rajan Saha Raju, Abdullah Al Nahid, Preonath Chondrow Dev,  Rashedul Islam. [VirusTaxo: Taxonomic classification of viruses from the genome sequence using k-mer enrichment
](https://www.sciencedirect.com/science/article/pii/S0888754322001598). Genomics, Volume 114, Issue 4, July 2022.
