#### 16S rRNA gene ####

# DECIPHER for 16S rRNA gene taxonomy assignment

library(DECIPHER)
library(tidyverse)

# SILVA r138 used as database
load("./SILVA_SSU_r138_2019.RData")

# load fasta file
fas <- "./clustered_mostabundant_16SrRNA.fas"
best_otus_seq <- readDNAStringSet(fas)

# taxonomy assignment
ids <- IdTaxa(best_otus_seq, trainingSet, strand = "both", processors = NULL, verbose = TRUE, threshold = 60)

# retrieve taxonomy
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks

# final taxonomy
tax <- as_tibble(taxid, rownames = "otu") %>% 
  separate("otu", 
           c("leaf","isolate"), # separate
           sep = "([\\|])", 
           remove = FALSE) %>%
  mutate(label = if_else(phylum == "Proteobacteria", true = class, false = phylum))

# select bacterial OTUs
bact_otus <- filter(tax, domain == "Bacteria") %>% 
  select(cluster)

# load OTU table
otus <- openxlsx::read_excel("./otu_table.xlsx")

# percentage of OTUs per samples
otus_bact_percent <- gather(otus, sample, count, 2:ncol(otus)) %>% 
  group_by(cluster) %>%
  mutate(cluster_seq_sum = sum(count)) %>% 
  filter(cluster_seq_sum > 1) %>% # remove global singletons
  ungroup %>%
  inner_join(bact_otus, by = "cluster") %>%  # remove non-bacterial OTUs
  group_by(sample) %>%
  mutate(seq_number_sample = sum(count)) %>% # count sum per sample
  ungroup %>%
  mutate(per = ((count / seq_number_sample) * 100))
 
#### ITS2 ####

# load OTU table
otus <- openxlsx::read_excel("./otu_table.xlsx")

# taxonomy from blastn
tax <- openxlsx::read_excel("./blastn.xlsx") %>% 
  separate("SEQ TITLE", 
           c("empty", "cluster", "mostabund_removed", "seq_count_removed"),
           sep = "([\\|])", 
           remove = FALSE) %>% 
  separate(Description, 
           c("lowest_tax", "accession", "sh", "refs", "taxonomy"), 
           sep = "([\\|])", 
           remove = FALSE) %>%
  separate(taxonomy, 
           c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
           sep = "([\\;])", 
           remove = FALSE) %>%
  select(-mostabund_removed, -seq_count_removed, -empty, -Accession) %>% 
  mutate(semifinal_tax = ifelse((Similarity > 97) & (Coverage > 90), lowest_tax, purrr::map_chr(genus, paste0, "_sp"))) %>%  # genus if sim or cov is low, otherwise species is assigned
  mutate(final_tax = ifelse((semifinal_tax == "g__unidentified_sp"), lowest_tax, semifinal_tax))

tax$final_tax=gsub("g__","", tax$final_tax)

# select fungal OTUs
fun_otus <- filter(tax, kingdom=="k__Fungi") %>% 
  select(cluster)

# percentage of OTUs per sample
otus_fun_percent <- gather(otus, sample, count, 2:ncol(otus)) %>% 
  group_by(cluster) %>%
  mutate(cluster_seq_sum = sum(count)) %>% 
  filter(cluster_seq_sum > 1) %>% # remove global singletons
  ungroup %>%
  inner_join(fun_otus, by = "cluster") %>%  # remove nonfungal OTUs
  group_by(sample) %>%
  dplyr::mutate(seq_number_sample = sum(count)) %>% # count sum per sample
  ungroup %>%
  dplyr::mutate(per = ((count / seq_number_sample) * 100))

