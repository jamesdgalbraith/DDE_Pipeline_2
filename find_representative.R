library(BSgenome)
library(plyranges)
library(tidyverse)
library(openxlsx)

# Get names of families
families <- getSheetNames("data/Supp-info_clades_CP_updt.xlsx")

# For each family
for(i in 1:length(families)){
  
  # Read in sheet of family
  sheet1 <- as_tibble(read.xlsx(xlsxFile = "data/Supp-info_clades_CP_updt.xlsx", sheet = families[i], colNames = F))
  # Name columns
  colnames(sheet1) <- c("clade", "seqnames")
  # select representatives
  sheet1 <- sheet1[(!is.na(sheet1$clade) & sheet1$clade != "Outgroup"),]
  
  # Read in family sequence
  in_seq <- readAAStringSet(paste0("data/", families[i], ".fasta"))
  # Select representaive sequences
  filtered_seq <- in_seq[names(in_seq) %in% sheet1$seqnames]
  # Add family name to sequence if not present
  names(filtered_seq) <- ifelse(grepl(families[i], names(filtered_seq)),
                                names(filtered_seq),
                                paste0(families[i], "_", names(filtered_seq)))
  
  # Write representatives to file
  writeXStringSet(filtered_seq, paste0("data/representative/", families[i], ".fasta"))

}
