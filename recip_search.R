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

# write to file as nucleotide and amino acid
writeXStringSet(tblastn_seq, paste0("out/plain_tblastn_initial_fastas/", "Glomeris_maerens", "_seq_nt.fasta"))
writeXStringSet(x = translate(tblastn_seq, if.fuzzy.codon = "solve"),
                filepath = paste0("out/plain_tblastn_initial_fastas/", "Glomeris_maerens", "_seq_aa.fasta"))

# search nucleotides against proteins
system(paste0("blastx -num_threads 8 -query out/plain_tblastn_initial_fastas/", "Glomeris_maerens",
              "_seq_nt.fasta -db data/representative/compiled_database/compiled.fasta -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue frames\" -out out/recip/", "Glomeris_maerens", "_nt.out"))

# read in recip tblastx
recip_blastx_nt <- read_tsv(paste0("out/recip/", "Glomeris_maerens", "_nt.out"),
                            col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue", "frames")
) %>%
  tidyr::separate(frames, into = c("qframe", "sframe"), sep = "/") 

# example of different reading frames
recip_blastx_nt_named %>%
  filter(qseqid =="55.93x_30964:7027-7457(-)", evalue <= 0.001) %>%
  group_by(seqnames) %>%
  arrange(seqnames, qstart)

# select best hits
recip_blastx_nt_named <- recip_blastx_nt %>%
  dplyr::filter(pident > 50, length >= 0.5*slen, evalue <= 0.001) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::slice(1) %>%
  filter(qframe != 1) %>%
  ungroup()  

tibble(classes = names(table(sub("_.*", "", recip_blastx_nt_named$seqnames))),
       n = as.integer(table(sub("_.*", "", recip_blastx_nt_named$seqnames))))


{system(paste0("blastp -num_threads 8 -query out/plain_tblastn_initial_fastas/", "Glomeris_maerens",
         "_seq_aa.fasta -db data/representative/compiled_database/compiled.fasta -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue\" -out out/recip/", "Glomeris_maerens", "_aa.out"))

recip_blastp_aa <- read_tsv(paste0("out/recip/", "Glomeris_maerens", "_aa.out"),
  col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue")
  )

recip_blastp_aa_named <- recip_blastp_aa %>%
  dplyr::filter(pident > 50, length >= 0.5*slen) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::slice(1) %>%
  ungroup()

tibble(classes = names(table(sub("_.*", "", recip_blastp_aa_named$seqnames))),
       n = as.integer(table(sub("_.*", "", recip_blastp_aa_named$seqnames))))
  
}