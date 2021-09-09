library(BSgenome)
library(plyranges)
library(tidyverse)

# read in blast output
tblastn_in <- read_tsv(paste0("out/tblastn/compiled_in_", "Glomeris_maerens", ".out"),
                       col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue"))

# read in genome
genome_seq <- readDNAStringSet(paste0("seq/", "Glomeris_maerens", ".fasta"))
names(genome_seq) <- gsub(" .*", "", names(genome_seq))

# filter blast output
tblastn_fixed <- tblastn_in %>%
  filter(length >= 0.5*qlen) %>%
  mutate(strand = ifelse(sstart < send, "+", "-"),
         start = ifelse(strand == "+", sstart, send),
         end = ifelse(strand == "+", send, sstart)) %>%
  as_granges()

# get filtered seq
tblastn_ranges <- GenomicRanges::reduce(tblastn_fixed)
tblastn_seq <- getSeq(genome_seq, tblastn_ranges)
names(tblastn_seq) <- paste0(seqnames(tblastn_ranges), ":", ranges(tblastn_ranges), "(", strand(tblastn_ranges), ")")

writeXStringSet(tblastn_seq, paste0("out/plain_tblastn_initial_fastas/", "Glomeris_maerens", "_seq_nt.fasta"))
writeXStringSet(x = translate(tblastn_seq, if.fuzzy.codon = "solve"),
                filepath = paste0("out/plain_tblastn_initial_fastas/", "Glomeris_maerens", "_seq_aa.fasta"))

system(paste0("blastp -num_threads 8 -query out/plain_tblastn_initial_fastas/", "Glomeris_maerens",
         "_seq_aa.fasta -db data/representative/compiled_database/compiled.fasta -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue\" -out out/recip/", "Glomeris_maerens", ".out"))

recip_blastp <- read_tsv(paste0("out/recip", "Glomeris_maerens", ".out"),
  col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue")
  )

recip_blastp %>%
  dplyr::group_by(qseqid)
