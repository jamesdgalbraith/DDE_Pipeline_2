#!/usr/bin/env Rscript

if(!"optparse" %in% installed.packages()){install.packages("optparse")}
if(!"tidyverse" %in% installed.packages()){install.packages("tidyverse")}
if(!"future.apply" %in% installed.packages()){install.packages("future.apply")}
if(!"BiocManager" %in% installed.packages()){install.packages("BiocManager")}
if(!"BSgenome" %in% installed.packages()){BiocManager::install("BSgenome")}
if(!"plyranges" %in% installed.packages()){BiocManager::install("plyranges")}

