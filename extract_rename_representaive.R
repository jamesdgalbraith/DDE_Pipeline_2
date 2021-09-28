library(BSgenome)
library(plyranges)
library(tidyverse)
library(openxlsx)
library(DECIPHER)

# Get names of families
families <- getSheetNames("data/Supp-info_clades_CP_updt.xlsx")

# For each family
for(i in 1:length(families)){
  
  # Read in sheet of family
  sheet1 <- as_tibble(read.xlsx(xlsxFile = "data/Supp-info_clades_CP_updt.xlsx", sheet = families[i], colNames = F))
  
  # Name columns
  colnames(sheet1) <- c("clade", "seqnames")
  
  # select representatives
  sheet1 <- sheet1[!is.na(sheet1$clade),]
  
  # Read in family sequence
  alignment <- readAAStringSet(paste0("data/alignments/og/", families[i], ".fasta"))

  # select representative sequence  
  alignment <- alignment[names(alignment) %in% sheet1$seqnames]
  
  # determine which is which
  names_alignment <- tibble(seqnames = names(alignment)) %>%
    inner_join(sheet1)
  tibble(names = names(table(sheet1$seqnames)), n = as.integer(table(sheet1$seqnames))) %>% arrange(-n)
  # names sequences accordingly
  names(alignment) <- paste0(families[i], "_", sub("Clade ", "", names_alignment$clade))
  
  # write unaligned to file
  writeXStringSet(x = RemoveGaps(alignment),
                  filepath = paste0("data/unaligned/", families[i], "_representative.fasta"))
  
  # Write representative aligned to file
  writeXStringSet(alignment, paste0("data/alignments//", families[i], ".fasta"))

}



