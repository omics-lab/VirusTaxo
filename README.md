### VirusTaxo

VirusTaxo provides taxonomic classification of virus sequences. VirusTaxo has an average accuracy of 93% at genus level across DNA and RNA viruses.

### Method limitation and interpretation

- VirusTaxo is trained on known virus sequences and designed to predict taxonomy of virus sequences. 

- Since VirusTaxo uses k-mer enrichment, non-virus sequences could be classified as virus due to random k-mer match. To use VirusTaxo, please make sure to remove non-virus sequences.    

- If your sample contains non-virus sequences, it is highly recommeded to filter our non-viral sequences using tools like [blast](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/find-data/sequence) or [DeepVirFinder](https://github.com/jessieren/DeepVirFinder). 

- To avoid contaimination with host sequences, please filter out those by mapping the reads to host reference genomes before using VirusTaxo. 

### Web application of VirusTaxo

- Web-based application of VirusTaxo is available at [Omics Lab](https://omics-lab.com/virustaxo) 


### Running VirusTaxo 
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

### Predict virus taxonomy from fasta file using prebuilt database

- Download prebuilt databse of VirusTaxo 
   - database.v2_2024 (recommended): download `vt_db_jan21_2024.tar.gz` from [here](https://drive.google.com/file/d/1gz0n5oHomWjpT0HXsrqh8hTLqmqiqgJs/view?usp=sharing).
   - database.v1_2022: download `vt_db_apr27_2022.tar.gz` from [here](https://drive.google.com/file/d/1j9rcFi6AMjA7tSqSizAQO7GpZw-brauZ/view?usp=sharing).
- Extract three database files using `tar –xvzf vt_db_jan21_2024.tar.gz`. 
   - DNA_RNA_18451_k20.pkl (recommended): database for both DNA and RNA viruses  
   - DNA_9384_k21.pkll : database for DNA viruses only
   - RNA_9067_k17.pkl : database for RNA viruses only

- Example of predicting virus taxonomy from the combined database 
   - Perform *de novo* assembly to generate `contig.fasta` file from your metavirome or metagenomic library
   - *De novo* assembly using [MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884) `megahit -1 file_R1.fq -2 file_R2.fq --min-contig-len 500 -o contig.fasta`
   - If needed, filter our non-viral sequences using tools like [blast](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/find-data/sequence) or [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) `dvf.py -i contig.fasta -o dvf.contig.fasta`  

```
python3 predict.py \
   --model_path /path/DNA_RNA_18451_k20.pkl \
   --seq ./contig.fasta
```

### Build your custom database

- Preparing a metadata file in `csv` format. The metadata file must contain two columns named `Id`  and `Genus`. For example:
```
Id,Genus
NC_004205.1,Seadornavirus
NC_038276.1,Ledantevirus
```

- Preparing a sequence file in `fasta` format. For example:

```
>NC_004205.1
GTAGAAATTTGTAAAGATTAACAATGTCGAGTTTAAAGGAACATAGGACTAATAAAGCGAATTCAAGAAA
CTTAATTCGTAGTCCAGATGAAGCTCCACCAACAGATAACAGTTTATTGAACAAAGGTGAAATATTAGCA
CTTACTTTCAGTGATGAATACATCAAATCAAAACTATTACTTGGTCCGAAACTGCAAGGTTTACCTCCTC
CATCACTTCCACCCAATTCGTACGGTTATCATTGCAATGGGTCGTTCGCCACCTATTTGCTAAGAGAATC
>NC_038276.1
CTTGAGAAACTTATTAGTCTATCAAGGTGGTTGTTTTTTCCCACATGAGTCAACAGACATCAAAATGTCT
GATAGAGTTCCCTTCCGTGTTGCTACCAAGCAACCTGTTAAGCCCATTCTTCCACAAGAGGAGACCCCAG
GACAATATCCAGCCGACTGGTTCAACACCCACAAAAATGAAAAGCCACGATTAGTCATTCCCTATAAAAT
CAAAGATATGGATTCTCTGAGAGGTATTGTTCGTGAAGGAATAGAAAAGGACAGCCTGGATGTTAAGGTG
```

 **N.B.** The ID of the sequences must present in metadata file.


 - DB Building cmd:
```
python3 build.py \
   --meta ./Dataset/RNA_meta.csv \
   --seq ./Dataset/RNA_seq.fasta \
   --k 17 \
   --saving_path ./model/RNA.pkl
```

 - Details of Parameters 
  
   - `meta`: Absolute or relative path of metadata file.
   - `seq`: Absolute or relative path of fasta sequence file.
   - `k` : It denotes the length of k-mer during database building.
   - `saving_path`: The program will save a pickle file (A DB File) in the mentioned path.

   
#### Predict virus taxonomy from fasta file using the custom database

- Prediction cmd with a sample input.fasta
```
 python3 predict.py \
   --model_path ./model/RNA.pkl \
   --seq ./Dataset/input.fasta
```

- Sample Output. Higher entropy (>0.5) is considered as `Unclassified` and lower entropy (≤0.5) provides the higher level of certainty at the genus level prediction.
- Genus name `Unknown` means the genus name is not assigned in the [ICTV classification](https://ictv.global/). 

```
Id      Length  Genus   Entropy
NC_004205.1     280     Seadornavirus   0.0
NC_038276.1     280     Ledantevirus    0.0
```

### Hierarchical classification 

[Find here the earlier version of VirusTaxo with hierarchical classification and the codes used in publication.](https://github.com/omics-lab/VirusTaxo_Hierarchical)

### Contact
Rashedul Islam, PhD (rashedul.gen@gmail.com)

### Citation

Rajan Saha Raju, Abdullah Al Nahid, Preonath Chondrow Dev,  Rashedul Islam. [VirusTaxo: Taxonomic classification of viruses from the genome sequence using k-mer enrichment
](https://www.sciencedirect.com/science/article/pii/S0888754322001598). Genomics, Volume 114, Issue 4, July 2022.
