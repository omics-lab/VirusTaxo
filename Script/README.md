## VirusTaxo: Taxonomic classification of viruses from metagenomic contigs

### 1. Running VirusTaxo 
#### Requirements 
- python >= 3.8
- Linux

#### Installation

 - Clone the repository

```
git clone https://github.com/omics-lab/VirusTaxo
```

 - Create Python Virtual Environment

```
cd VirusTaxo
python3 -m venv environment
source ./environment/bin/activate
```
 - Install Python Packages

```
pip install -r requirements.txt
```

### 2. Predict virus taxonomy using prebuilt database

#### Step-1: Download the [latest]() prebuilt databse of VirusTaxo 

```
gdown "https://drive.google.com/uc?id=1gz0n5oHomWjpT0HXsrqh8hTLqmqiqgJs"

# Extract db files
tar -xvzf vt_db_jan21_2024.tar.gz
```

-  Database files:

| Database file | Count | Description |
|----------|----------|----------|
| Family_database.pkl  | 242 Family  | k-mer database for Family level prediction   |
| Genus_database.pkl  | 1933 Genus  | k-mer database for Genus level prediction   |
| Species_database.pkl  | 8528 Species  | k-mer database for Species level prediction   |
| sequences.fasta  | 12613 genome | Complete genome sequences used to build database |
| metadata.csv  | 12613 accession | Metadata associated with the dataset used to build database |

#### Step-2: Assemble the metagenomic contigs from your metavirome or metagenomic library 
   - Perform *de novo* assembly using [MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884): 

```
megahit -1 file_R1.fq -2 file_R2.fq --min-contig-len 500 -o contig.fasta
```

#### Step-3: Predict taxonomy using `predict.py`

- Usage

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

### 3. Interpretation of output

- Example output 

```
Accession    Query Seq Length  Family           Family Entropy  Family Enrichment  Genus           Genus Entropy  Genus Enrichment  Species                   Species Entropy  Species Enrichment  Valid
AC 000001.1  33034             Adenoviridae     0.0             0.786              Mastadenovirus  0.0            0.782             Ovine mastadenovirus A    0.0              0.745               Yes
AC 000002.1  34446             Adenoviridae     0.0             0.845              Mastadenovirus  0.0            0.842             Bovine mastadenovirus B   0.0              0.811               Yes
AC 000011.1  36519             Adenoviridae     -0.0            0.57               Mastadenovirus  -0.0           0.563             Human mastadenovirus E    -0.0             0.353               Yes
AC 000189.1  34094             Adenoviridae     -0.0            0.81               Mastadenovirus  -0.0           0.803             Porcine mastadenovirus A  -0.0             0.742               Yes
NC 000852.5  330611            Phycodnaviridae  -0.0            0.121              Chlorovirus     -0.0           0.115             Unclassified              -0.0             0.014               Yes
NC 000855.1  11158             Unclassified     0.0             0.011              Unclassified    0.0            0.006             Unclassified              0.09             0.002               Yes
NC 000867.1  10079             Unclassified     -0.0            0.018              Unclassified    -0.0           0.01              Unclassified              0.089            0.002               Yes
NC 000899.1  45063             Adenoviridae     0.0             0.85               Aviadenovirus   0.0            0.844             Fowl aviadenovirus D      0.0              0.772               Yes
NC 000939.2  4415              Tombusviridae    -0.0            0.061              Aureusvirus     -0.0           0.055             Aureusvirus dioscoreae    -0.0             0.053               Yes
```

- In the taxonomic rank column 

   - `Unclassified`: `Entropy` (default >= 0.5) or `Enrichment` (default <= 0.05 for Family and Genus; default <= 0.80 for Species classification) is outside of cutoff. 

   - Lower `Entropy` (such as ≤=0.5) provides the higher level of prediction certainty. You can decrease `Entropy` cutoff for better prediction. 

   - Higher `Enrichment_Score` (such as >= 0.8) provides the higher level of prediction certainty. You can increase `Enrichment_Score` cutoff for better prediction. `Enrichment_Score` is the total number of k-mers mapped to the genera divided by total number of k-mers in the query sequence.

   - `Valid` column is `Yes` if the prediction follows known taxonomic ranks, otherwise `No`. 

### 4. Prediction accurary of VirusTaxo
To check accuracy, 12,613 complete virus genomes were used. In 5-fold cross-validation, 80% of the sequences were randomly chosen to create the database, and the other 20% were used to calculate the accuracy shown in the table below:

|  | Family | Genus |  Species | 
|----------|----------|----------|----------|  
| Accurate prediction  | 97%  | 97%   | 87%   | 
| Wrong prediction  | 3%  | 3% | 13% | 
| Unclassified  | 45%  | 45%  | 55%  | 
| k-mer size (bp) | 15 | 15 | 20 | 
| Enrichment cutoff  | >=0.05 | >=0.05 | >=0.80 | 
| Entropy cutoff  | <=0.50 | <=0.50 | <=0.50 | 

### 5. Build your custom database

- Preparing a metadata file in `csv` format. The metadata file must contain columns named `Accession`, `Family`, `Genus` and `Species`. Example of metadata file is [here](./Dataset/RNA_meta.csv):

- The sequnce `Accession` must match with the metadata `Accession`. Example of input fasta file is [here](./Dataset/RNA_seq.fasta).

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


### 6. Method limitation and interpretation

- VirusTaxo is trained on known virus sequences and designed to predict taxonomy of virus sequences. 

- Since VirusTaxo uses k-mer enrichment, non-virus sequences could be classified as virus due to random k-mer match. To use VirusTaxo, please make sure to remove non-virus sequences.    

- To avoid contaimination with host sequences, please filter out those by mapping the reads to host reference genomes before using VirusTaxo. 

- If your sample contains non-virus sequences, it is recommeded to filter out non-viral sequences using tools like [blast](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/find-data/sequence) or [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) `dvf.py -i contig.fasta -o ./`. 

### 7. Version history 

| VirusTaxo | Database  | Data Date     | Sequences | Download |
|----------|----------|----------|----------|----------|
| [v2]() Family, Genus, Species prediction | database.v3_2024  | Jan21_2024  | DNA=9384 &  RNA=9067  | [here]()  |
| [v1](Script/v1) Genus prediction | database.v2_2024  | Jan21_2024  | DNA=9384 &  RNA=9067  | [here](https://drive.google.com/file/d/1gz0n5oHomWjpT0HXsrqh8hTLqmqiqgJs/view?usp=sharing)  |
| [v1](Script/v1) Genus prediction | database.v1_2022  | Apr27_2022  | DNA=4421 &  RNA=2529  | [here](https://drive.google.com/file/d/1j9rcFi6AMjA7tSqSizAQO7GpZw-brauZ/view?usp=sharing)  |
| [Used in manuscript](https://github.com/omics-lab/VirusTaxo_Hierarchical) | database.v1_2022  | Apr27_2022  | DNA=4421 &  RNA=2529  | [here](https://drive.google.com/file/d/1j9rcFi6AMjA7tSqSizAQO7GpZw-brauZ/view?usp=sharing) |

### 8. Contact
Rashedul Islam, PhD (rashedul.gen@gmail.com)

### 9. Citation

Rajan Saha Raju, Abdullah Al Nahid, Preonath Chondrow Dev,  Rashedul Islam. [VirusTaxo: Taxonomic classification of viruses from the genome sequence using k-mer enrichment
](https://www.sciencedirect.com/science/article/pii/S0888754322001598). Genomics, Volume 114, Issue 4, July 2022.
