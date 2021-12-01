#!/usr/bin/env Rscript

if(!"optparse" %in% installed.packages()){install.packages("optparse")}
if(!"tidyverse" %in% installed.packages()){install.packages("tidyverse")}
if(!"future.apply" %in% installed.packages()){install.packages("future.apply")}
if(!"openxlsx" %in% installed.packages()){install.packages("openxlsx")}
if(!"BiocManager" %in% installed.packages()){install.packages("BiocManager")}
if(!"BSgenome" %in% installed.packages()){BiocManager::install("BSgenome")}
if(!"ORFik" %in% installed.packages()){BiocManager::install("ORFik")}
if(!"plyranges" %in% installed.packages()){BiocManager::install("plyranges")}
if(!"DECIPHER" %in% installed.packages()){BiocManager::install("DECIPHER")}

library(tidyverse)

fastas <- read_tsv("fastas.txt", col_names = "fasta")

completed <- read_tsv("completed_genomes.txt", col_names = c("fasta", "species"))

completed[!completed$fasta %in% fastas$fasta,]

