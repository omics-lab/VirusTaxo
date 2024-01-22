# data source https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SourceDB_s=RefSeq
# download: select columns for the table and use downlowd button at the top to get seq and meta
# date: 21 Jan 2024

library(tidyverse)


meta = read.csv("./Dataset/sequences_20240122_8132761_with_version.csv")
head(meta)
colnames(meta)

meta2 = meta %>% filter(Nuc_Completeness == "complete") 

meta2$molecule = ifelse(grepl("DNA", meta2$Molecule_type), "DNA",
                        ifelse(grepl("RNA", meta2$Molecule_type), "RNA", "other"))

meta3 = meta2 %>% filter(molecule != "other")
colnames(meta3)
meta4 = meta3 %>% select("Accession", "molecule", "Molecule_type", "Organism_Name", "Species", "Genus", "Family", "Length")

# replace missing cells
meta5 = meta4
meta5[meta5 == ""] <- "Unknown"

# replace colname
colnames(meta5)[1] <- "Id"

dna_rna = meta5
dna = meta5 %>% filter(molecule == "DNA")
rna = meta5 %>% filter(molecule == "RNA")

write.csv(dna_rna, "DNA_RNA_18451_meta.csv")
write.csv(dna, "DNA_9384_meta.csv")
write.csv(rna, "RNA_9067_meta.csv")

# end




  
  
  
  
  