#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})


# parse input variables
option_list = list(
  make_option(c("-s", "--species_name"), type="character", default=NULL, 
              help="species name", metavar="character"),
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.na(opt$species_name)) {
  stop("Species name is needed")
} else {
  # set species names
  species_name <- opt$species_name
  message(species_name)
}

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

# read in recip tblastx
recip_blast <- read_tsv(paste0("out/recip/", species_name, ".out"),
                            col_names = c("seqnames", "sseqid", "pident", "length", "qstart", "qend", "qlen",
                                          "sstart", "send", "slen")
)

# read in sequence
both_seq <- readAAStringSet(filepath = paste0("out/plain_tblastn_initial_fastas/", species_name, "_seq.fasta"))

# filter based on length
recip_blast_filtered <- recip_blast %>%
  group_by(seqnames) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(sseqid = sub("_.*", "", sseqid))  %>%
  filter(length >= 0.8*slen, length <= 1.2*slen, qlen/slen<=1.25, qlen/slen>=0.75)

# make tibble to name sequences
classified_tbl <- tibble(seqnames = names(both_seq)) %>%
  inner_join(recip_blast_filtered)

# select classified sequences and name
classified_seq <- both_seq[names(both_seq) %in% classified_tbl$seqnames]
names(classified_seq) <- paste0(classified_tbl$seqnames, "#", classified_tbl$sseqid)

# create outdir if missing
if(!dir.exists("out/classified_tnps/")){dir.create("out/classified_tnps/")}

# write compiled
writeXStringSet(classified_seq, paste0("out/classified_tnps/", species_name, "_compiled.fasta"))

# write individual families
families <- base::unique(classified_tbl$sseqid)

for(i in 1:length(families)){
  
  # create outdir if missing
  if(!dir.exists(paste0("out/classified_tnps/", families[i]))){dir.create(paste0("out/classified_tnps/", families[i]))}
  
  # write to file
  writeXStringSet(x = classified_seq[classified_tbl$sseqid == families[i]],
                  filepath = paste0("out/classified_tnps/", families[i], "/", species_name, "_", families[i], ".fasta"))

}
