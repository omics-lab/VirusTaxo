## VirusTaxo: Taxonomic classification of viruses from metagenomic contigs

VirusTaxo provides taxonomic classification of virus sequences. VirusTaxo has an average accuracy of 93% at genus level across DNA and RNA viruses.

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

- Download prebuilt databse of VirusTaxo `database.v2_2024` from [here](https://drive.google.com/file/d/1gz0n5oHomWjpT0HXsrqh8hTLqmqiqgJs/view?usp=sharing).
- Extract three database files using `tar –xvzf vt_db_jan21_2024.tar.gz`. 


| db file | Molecule | Usage |
|----------|----------|----------|
| DNA_RNA_18451_k20.pkl  | DNA & RNA  | Samples containing both DNA & RNA viruses   |
| DNA_9384_k21.pkl  | DNA  | Samples containing DNA viruses only |
| RNA_9067_k17.pkl  | RNA  | Samples containing RNA viruses only |


- Assemble the metagenomic contigs from your metavirome or metagenomic library. 
   - Perform *de novo* assembly using [MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884): 

   ```megahit -1 file_R1.fq -2 file_R2.fq --min-contig-len 500 -o contig.fasta```

- Run with the sample `contig.fasta` file

```
python3 predict.py \
   --model_path /path/DNA_RNA_18451_k20.pkl \
   --seq ./Dataset/contig.fasta
```

- Example output for 4 query sequences

```
Id              Length  Genus           Entropy Enrichment_Score
QuerySeq-1      219     NoHit           1.0     0
QuerySeq-2      720     Betacoronavirus 0.0     0.974
QuerySeq-3      1540    Unknown_genus   0.527   0.0020
QuerySeq-4      1330    Lentivirus      0.0     0.991
```

- Example output after filtering the query sequences by `Entropy` <=0.5 and `Enrichment_Score` >=0.8 (**Recommended**)

```
Id              Length  Genus           Entropy Enrichment_Score
QuerySeq-2      720     Betacoronavirus 0.0     0.974
QuerySeq-4      1330    Lentivirus      0.0     0.991
```

### 3. Interpretation of output
- `NoHit` means there's no k-mer overlap between the query and database.
- Lower `Entropy` (such as ≤=0.5) provides the higher level of prediction certainty. You can decrease `Entropy` cutoff for better prediction. 
   - We recommend to filter out the query sequences with `Entropy` cutoff of ≤0.5. 
- Higher `Enrichment_Score` (such as >= 0.8) provides the higher level of prediction certainty. You can increase `Enrichment_Score` cutoff for better prediction. `Enrichment_Score` is the total number of k-mers mapped to the genera divided by total number of k-mers in the query sequence.
   - We recommend to filter out the query sequences with `Enrichment_Score` cutoff of >=0.8. 
- Genus name `Unknown` means the genus name is not assigned in the [ICTV classification](https://ictv.global/). 


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
   - `saving_path`: Path to save a pickle file (A DB File).


### 5. Method limitation and interpretation

- VirusTaxo is trained on known virus sequences and designed to predict taxonomy of virus sequences. 

- Since VirusTaxo uses k-mer enrichment, non-virus sequences could be classified as virus due to random k-mer match. To use VirusTaxo, please make sure to remove non-virus sequences.    

- If your sample contains non-virus sequences, it is highly recommeded to filter out non-viral sequences using tools like [blast](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/find-data/sequence) or [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) `dvf.py -i contig.fasta -o ./`. 

- To avoid contaimination with host sequences, please filter out those by mapping the reads to host reference genomes before using VirusTaxo. 


### 6. Database versions

| Version  | Date     | Sequences | Download |
|----------|----------|----------|----------|
| database.v2_2024  | Jan21_2024  | DNA=9384 &  RNA=9067  | [here](https://drive.google.com/file/d/1gz0n5oHomWjpT0HXsrqh8hTLqmqiqgJs/view?usp=sharing)  |
| database.v1_2022  | Apr27_2022  | DNA=4421 &  RNA=2529  | [here](https://drive.google.com/file/d/1j9rcFi6AMjA7tSqSizAQO7GpZw-brauZ/view?usp=sharing)  |

### 7. Hierarchical classification 

[Find here the earlier version of VirusTaxo with hierarchical classification and the codes used in publication.](https://github.com/omics-lab/VirusTaxo_Hierarchical)

### 8. Contact
Rashedul Islam, PhD (rashedul.gen@gmail.com)

### 9. Citation

Rajan Saha Raju, Abdullah Al Nahid, Preonath Chondrow Dev,  Rashedul Islam. [VirusTaxo: Taxonomic classification of viruses from the genome sequence using k-mer enrichment
](https://www.sciencedirect.com/science/article/pii/S0888754322001598). Genomics, Volume 114, Issue 4, July 2022.
