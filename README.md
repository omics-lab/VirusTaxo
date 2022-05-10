### VirusTaxo

For taxonomic classification of viruses from metagenomic sequences, VirusTaxo builds database from diverse (e.g., 402 DNA and 280 RNA) genera of viruses. VirusTaxo has an average accuracy of 93% at genus level prediction across DNA and RNA viruses.

### Web application of virus taxo

- Web-based application of VirusTaxo is available at [Omics Lab](https://omics-lab.com/virustaxo) 


### Running VirusTaxo from the command line
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

#### Predict virus taxonomy from fasta file using prebuilt database

- Download prebuilt databse of VirusTaxo `vt_db_apr27_2022.tar.gz` from [here](https://drive.google.com/file/d/1j9rcFi6AMjA7tSqSizAQO7GpZw-brauZ/view?usp=sharing).
- Extract three database files using `tar –xvzf vt_db_apr27_2022.tar.gz`. 
   - vt_db_all_virus_kmer_20.pkl  (combined database for DNA and RNA viruses)  
   - vt_db_dna_virus_kmer_21.pkl  (database for DNA viruses)
   - vt_db_rna_virus_kmer_17.pkl  (database for RNA viruses)

- Example of predicting virus taxonomy from the combined database 
   - Perform de novo assembly to generate `input_contig.fasta` file from your metagenomic library
   - Example for de novo assembly using [MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884) `megahit -1 pe_1.fq -2 pe_2.fq -o out`

```
python3 predict.py \
   --model_path /path/vt_db_all_virus_kmer_20.pkl \
   --seq ./input_contig.fasta
```

#### Build custom database

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

Sample input.fasta

```
>k141_107617
CCATTAGTACATTTCTGTGCATACAATATGATAGATCTGATCAAGATCATTATGCTCAAT
ATCAATTAGTCAATCAAGTCATGTGTAATTATAACTACTCTTAATATTGCACTTGTACTC
ATTGTCCTGACTGGCAACTAAAAACATACATCTGCTTCTTAGACTATGTATTTACATGTA
CATGGTTGTGAATTGAATATGAAAGTGCTACAGTATCATAGTAACATGTACCATATTTAT
CAGTGCAACTAAATTCTACATGAGATAGGGTAATATATTTATATTTGTTACATGTATTTT
ACCTTGTAC
>k141_21524
CAGGTGCAGATACATATAGTCTAAGCTCCTCTTGATTAGCAATTACAGCAGCAGGAATAG
CAGCATATACTAAAGCTAATTTAGCAAGTACATTACTAGCATTAACAGCTACAGGAGATG
CGATATCAATTACATTAGCAGAATCAGCTACTAAAGATACTTTGTATCCATCACATAAAG
CAAGAGCAGGAGTACCTGAATTGACATCACCTGACCAACGCAATTTCTCTACATTCTC
```

- Prediction cmd
```
 python3 predict.py \
   --model_path ./model/RNA.pkl \
   --seq ./input.fasta
```

- Sample Output

```
Id              Length  Genus           Entropy
k141_107617     309     Unclassified    1.0
k141_21524      238     Tobamovirus     0.0
```

   
### Citation

Rajan Saha Raju, Abdullah Al Nahid, Preonath Shuvo,  Rashedul Islam. [VirusTaxo: Taxonomic classification of virus genome using multi-class hierarchical classification by k-mer enrichment](https://www.biorxiv.org/content/10.1101/2021.04.29.442004v1.full). bioRxiv, April 30, 2021.
