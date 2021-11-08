#!/usr/bin/Rscript

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

if(!dir.exists("out/tree_data/")){dir.create("out/tree_data/")}

## AMINO ACID ##
# read in families found
families_found <- read_tsv("out/orfs/aa/families_found.txt", col_names = "family")

i=7

# read in family of interest
family_clustered <- read_tsv(paste0("out/orfs/aa/", families_found$family[i], ".fasta.cd.clstr.tsv")) %>%
  mutate(species = sub("#.*", "", id),
         clstr_cov = as.double(sub("%", "", clstr_cov)),
         clstr_iden = as.double(sub("%", "", clstr_iden)))

# identify clusters with more than 1 member
larger_clusters <- as_tibble(as.data.frame(table(family_clustered$clstr))) %>%
  mutate(Var1 = as.integer(Var1)) %>%
  filter(Freq > 1)

# select clusters with more than 1 member
larger_clusters <- family_clustered[family_clustered$clstr %in% larger_clusters$Var1,]

# identify clusters containing multiple species
uniq_species_clstr <- base::unique(tibble(clstr = larger_clusters$clstr, species = larger_clusters$species))
multiple_species_clstrs <- tibble(clstr = names(table(uniq_species_clstr$clstr)), Freq = as.integer(table(uniq_species_clstr$clstr))) %>%
  filter(Freq > 2) %>%
  mutate(clstr = as.integer(clstr)) %>%
  arrange(-Freq)

# include number of species with cluster info
multiple_species_clstrs <- larger_clusters %>%
  inner_join(multiple_species_clstrs)

# select potential genes based on mean clstr identity and number of species
potential_genes <- multiple_species_clstrs[multiple_species_clstrs$clstr_size <= 1.2*multiple_species_clstrs$Freq,] %>%
  arrange(-Freq, clstr_size, clstr, species) %>%
  filter(Freq > 3) %>%
  group_by(clstr) %>%
  mutate(mean_iden = base::mean(clstr_iden)) %>%
  ungroup()

View(potential_genes)

# vector of clusters which are potentiallly genes
potential_gene_clstrs <- base::unique(potential_genes$clstr)

larger_clusters_reps <- larger_clusters %>%
  filter(clstr_rep == 1)

larger_cluster_aa_seq <- readAAStringSet(paste0("out/orfs/aa/", families_found$family[i], ".fasta"))

larger_cluster_aa_seq <- larger_cluster_aa_seq[names(larger_cluster_aa_seq) %in% larger_clusters_reps$id]

larger_cluster_aa_names <- tibble(id = names(larger_cluster_aa_seq)) %>%
  inner_join(larger_clusters_reps) %>%
  mutate(names = paste0(families_found$family[i], "_cluster_", clstr))

names(larger_cluster_aa_seq) <- paste0(larger_cluster_aa_names$names)

# writeXStringSet(x = larger_cluster_aa_seq, filepath = paste0(families_found$family[i]))


# ## NUCLEOTIDE ##
# 
# # read in families found
# families_found <- read_tsv("out/orfs/nt/families_found.txt", col_names = "family")
# 
# i=7
# 
# # read in family of interest
# family_clustered <- read_tsv(paste0("out/orfs/nt/", families_found$family[i], ".fasta.cd.clstr.tsv")) %>%
#   mutate(species = sub("#.*", "", id),
#          clstr_cov = as.double(sub("%", "", clstr_cov)),
#          clstr_iden = as.double(sub("%", "", clstr_iden)))
# 
# # identify clusters with more than 1 member
# larger_clusters <- as_tibble(as.data.frame(table(family_clustered$clstr))) %>%
#   mutate(Var1 = as.integer(Var1)) %>%
#   filter(Freq > 1)
# 
# # select clusters with more than 1 member
# larger_clusters <- family_clustered[family_clustered$clstr %in% larger_clusters$Var1,]
# 
# # identify clusters containing multiple species
# uniq_species_clstr <- base::unique(tibble(clstr = larger_clusters$clstr, species = larger_clusters$species))
# multiple_species_clstrs <- tibble(clstr = names(table(uniq_species_clstr$clstr)), Freq = as.integer(table(uniq_species_clstr$clstr))) %>%
#   filter(Freq > 1) %>%
#   mutate(clstr = as.integer(clstr))
# 
# # include number of species with cluster info
# multiple_species_clstrs <- larger_clusters %>%
#   inner_join(multiple_species_clstrs)
# 
# # select potential genes based on mean clstr identity and number of species
# potential_genes <- multiple_species_clstrs[multiple_species_clstrs$clstr_size <= 1.5*multiple_species_clstrs$Freq,] %>%
#   filter(Freq > 2) %>%
#   arrange(-Freq, clstr_size, clstr, species) %>%
#   group_by(clstr) %>%
#   mutate(mean_iden = base::mean(clstr_iden)) %>%
#   ungroup() %>%
#   filter(mean_iden > 85)
# 
# # vector of clusters which are potentiallly genes
# potential_gene_clstrs <- base::unique(potential_genes$clstr)
# 
# larger_clusters_reps <- larger_clusters %>%
#   filter(clstr_rep == 1)
# 
# larger_cluster_nt_seq <- readAAStringSet(paste0("out/orfs/nt/", families_found$family[i], ".fasta"))
# 
# larger_cluster_nt_seq <- larger_cluster_nt_seq[names(larger_cluster_nt_seq) %in% larger_clusters_reps$id]
# 
# larger_cluster_nt_names <- tibble(id = names(larger_cluster_nt_seq)) %>%
#   inner_join(larger_clusters_reps) %>%
#   mutate(names = paste0(families_found$family[i], "_cluster_", clstr))
# 
# names(larger_cluster_nt_seq) <- paste0(larger_cluster_nt_names$names)
# 
# writeXStringSet(x = larger_cluster_nt_seq, filepath = paste0(families_found$family[i]))
