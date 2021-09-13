library(BSgenome)
library(plyranges)
library(tidyverse)

# read in recip tblastx
recip_blast <- read_tsv(paste0("out/recip/", "Glomeris_maerens", ".out"),
                            col_names = c("seqnames", "qseqid", "pident", "length", "qstart", "qend", "qlen",
                                          "sstart", "send", "slen", "evalue", "frames")
) %>%
  tidyr::separate(frames, into = c("qframe", "sframe"), sep = "/") 

both_seq <- readAAStringSet(filepath = paste0("out/plain_tblastn_initial_fastas/", genome_name, "_seq.fasta"))

recip_blast_sliced <- recip_blast %>%
  group_by(seqnames) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(qseqid = sub("_.*", "", qseqid))

recip_blast_filtered <- recip_blast %>%
  group_by(seqnames) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(qseqid = sub("_.*", "", qseqid))  %>%
  # filter(pident >= 50) %>%
  filter(length >= 0.8*slen, length <= 1.2*slen)

table(recip_blast_filtered$qseqid)

classified_tbl <- tibble(seqnames = names(both_seq)) %>%
  inner_join(recip_blast_filtered)

classified_seq <- both_seq[names(both_seq) %in% classified_tbl$seqnames]
names(classified_seq) <- paste0(classified_tbl$seqnames, "#", classified_tbl$qseqid)

family_oi <- "Tc1marPlm"

of_interest <- classified_seq[classified_tbl$qseqid == family_oi]
writeXStringSet(of_interest, "temp/temp.fasta")

ggplot(recip_blast_filtered, aes(x = pident, y = length, colour = qseqid)) + geom_point()
