#!/usr/bin/env Rscript

library(tidyverse)
library(BSgenome)
library(plyranges)
library(openxlsx)

families <- getSheetNames("data/Supp-info_clades_CP_updt.xlsx")

for (i in 1:length(families)) {

  sheet1 <- as_tibble(read.xlsx(xlsxFile = "data/Supp-info_clades_CP_updt.xlsx", sheet = families[i], colNames = F)) %>%
    fill(X1)
  
  colnames(sheet1) <- c("clade", "seqnames")
  
  family_seq <- readAAStringSet(paste0("data/unaligned/", families[i], ".fasta"))
  
  family_info <- tibble(seqnames = names(family_seq), width = width(family_seq))
  
  longest_clades <- inner_join(sheet1, family_info) %>%
    group_by(clade) %>%
    arrange(-width) %>%
    slice(1) %>%
    ungroup()
  
  ggplot(longest_clades, aes(width)) + geom_density() + scale_x_continuous()
  
  longest_clades %>%
    arrange(width)
  
  writeXStringSet(family_seq[names(family_seq) %in% longest_clades$seqnames], paste0("data/unaligned/", families[i], "_representative.fasta"))
  
}
