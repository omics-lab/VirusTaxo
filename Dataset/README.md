### VirusTaxo Dataset

#### data download
- [NCBI data](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SourceDB_s=RefSeq) Select columns of the table and use downlowd button at the top to get sequence and metadata
- [metadata](./sequences_20240122_8132761_with_version.csv) and [fasta sequence](https://drive.google.com/file/d/1FBPgDFVvkIrfJ6XwEbpCPZesiaAbF67d/view?usp=sharing) are available to download
- download date: 21 Jan 2024
- All covid19 assembly by [Islam R et al, 2021](https://academic.oup.com/bib/article/22/5/bbab102/6210065) are available [here](https://drive.google.com/file/d/1fEPE3dcSMU4Ojq4T_C67owoIGfLwky6Y/view?usp=sharing)

```
.
├── asm_islam.et.al_6648_covid.fasta (covid19 assembly)
├── DNA_9384_meta.csv (DNA metadata latest)
├── DNA_RNA_18451_meta.csv (DNA and RNA metadata latest)
├── README.md 
├── RNA_9067_meta.csv (RNA metadata latest)
├── RNA_meta.csv (RNA metadata previous)
├── RNA_seq.fasta (RNA virus sequence previous)
└── RNA_0.5k_test.pkl (test database with 500 RNA sequence)
```

### Contact
Rashedul Islam, PhD (rashedul.gen@gmail.com)

### Citation

Rajan Saha Raju, Abdullah Al Nahid, Preonath Chondrow Dev,  Rashedul Islam. [VirusTaxo: Taxonomic classification of viruses from the genome sequence using k-mer enrichment
](https://www.sciencedirect.com/science/article/pii/S0888754322001598). Genomics, Volume 114, Issue 4, July 2022.
