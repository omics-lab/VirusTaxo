### Requirements
 - python >= 3.8
 - Linux

### Installation
 - Cloning the repository
```
git clone https://github.com/omics-lab/VirusTaxo_v1
```
 - Creation of Python Virtual Environment
```
cd VirusTaxo_v1
python3 -m venv environment
source ./environment/bin/activate
```
 - Installation of Python Packages
```
pip install -r requirements.txt
```

### Build Custom Database

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

   
### Predict from fasta file using Prebuilt DB

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