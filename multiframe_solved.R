library(BSgenome)
library(plyranges)
library(tidyverse)
library(ips)

# read in blast output
tblastn_in <- read_tsv(paste0("out/tblastn_pssm/compiled_in_", "Glomeris_maerens", ".out"),
                       col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue", "frames"))

# read in genome
genome_seq <- readDNAStringSet(paste0("seq/", "Glomeris_maerens", ".fasta"))
names(genome_seq) <- gsub(" .*", "", names(genome_seq))

# filter blast output
tblastn_fixed <- tblastn_in %>%
  dplyr::filter(length >= 0.5*qlen) %>%
  tidyr::separate(frames, into = c("qframe", "sframe"), sep = "/") %>%
  dplyr::mutate(strand = ifelse(sstart < send, "+", "-"),
         start = ifelse(strand == "+", sstart, send),
         end = ifelse(strand == "+", send, sstart),
         qframe = as.integer(qframe), sframe = as.integer(sframe))
  
# Split into separate into each frame
tblastn_fwd <- tblastn_fixed %>%
  dplyr::filter(strand == "+") %>%
  plyranges::as_granges()

# Seperate into individual frames
tblastn_fwd_1 <- GenomicRanges::reduce(tblastn_fwd[tblastn_fwd$sframe == 1]) %>% dplyr::mutate(sframe = 1)
tblastn_fwd_2 <- GenomicRanges::reduce(tblastn_fwd[tblastn_fwd$sframe == 2]) %>% dplyr::mutate(sframe = 2)
tblastn_fwd_3 <- GenomicRanges::reduce(tblastn_fwd[tblastn_fwd$sframe == 3]) %>% dplyr::mutate(sframe = 3)

# Reduce/merge (will contain different frames)
tblastn_fwd_reduced <- GenomicRanges::reduce(tblastn_fwd)

# Overlap frames and name columns
tblastn_fwd_combined <- rbind(plyranges::pair_overlaps(tblastn_fwd_reduced, tblastn_fwd_1),
      plyranges::pair_overlaps(tblastn_fwd_reduced, tblastn_fwd_2)) %>%
  rbind(plyranges::pair_overlaps(tblastn_fwd_reduced, tblastn_fwd_3))
colnames(tblastn_fwd_combined) <- c("joined", "framed", "frame")

# Determine those in single frames and make ranges
single_frames_fwd <- tibble(joined = names(table(tblastn_fwd_combined$joined)), n = as.integer(table(tblastn_fwd_combined$joined))) %>% filter(n==1) %>%
  tidyr::separate(joined, into = c("seqnames", "ranges", "strand"), sep = ":") %>%
  tidyr::separate(ranges, into = c("start", "end"), sep = "-") %>%
  dplyr::mutate(start = as.integer(start), end = as.integer(end)) %>%
  plyranges::as_granges()

# Determine those in multiple frames and make ranges
multiple_frames_fwd <- tibble(joined = names(table(tblastn_fwd_combined$joined)), n = as.integer(table(tblastn_fwd_combined$joined))) %>%
  filter(n==2) %>%
  tidyr::separate(joined, into = c("seqnames", "ranges", "strand"), sep = ":") %>%
  tidyr::separate(ranges, into = c("start", "end"), sep = "-") %>%
  dplyr::mutate(start = as.integer(start), end = as.integer(end)) %>%
  plyranges::as_granges()

# Select those in each frame and determine match
tblastn_fwd_1a <- pair_overlaps(tblastn_fwd_1,multiple_frames_fwd)$granges.x %>% 
  dplyr::mutate(joined = as.character(pair_overlaps(tblastn_fwd_1,multiple_frames_fwd)$granges.y), frame = 1)

tblastn_fwd_2a <- pair_overlaps(tblastn_fwd_2,multiple_frames_fwd)$granges.x %>% 
  dplyr::mutate(joined = as.character(pair_overlaps(tblastn_fwd_2,multiple_frames_fwd)$granges.y), frame = 2)

tblastn_fwd_3a <- pair_overlaps(tblastn_fwd_3,multiple_frames_fwd)$granges.x %>% 
  dplyr::mutate(joined = as.character(pair_overlaps(tblastn_fwd_3,multiple_frames_fwd)$granges.y), frame = 3)

tblastn_fwd_1_2a <- c(tblastn_fwd_1a[tblastn_fwd_1a$joined %in% tblastn_fwd_2a$joined],
                      tblastn_fwd_2a[tblastn_fwd_2a$joined %in% tblastn_fwd_1a$joined])

list_to_loop <- tblastn_fwd_1a[tblastn_fwd_1a$joined %in% tblastn_fwd_2a$joined]$joined
for(i in 1:list_to_loop){
  ranges_a <- tblastn_fwd_1_2a %>%
    filter(joined == list_to_loop[i])
  
  seq_a <- getSeq(genome_seq, ranges_a)
  
  
}






# Determine those in all three frames and make ranges
tblastn_fwd_1_2 <- pair_overlaps(tblastn_fwd_1, tblastn_fwd_2)
tblastn_fwd_1_2_3 <- pair_overlaps(GenomicRanges::reduce(c(tblastn_fwd_1_2$granges.x, tblastn_fwd_1_2$granges.y)), tblastn_fwd_3)
if(nrow(tblastn_fwd_1_2_3)>0){
  tblastn_fwd_1_2_3 <- GenomicRanges::reduce(c(tblastn_fwd_1_2_3$granges.x, tblastn_fwd_1_2_3$granges.y))
  }

findOverlaps(tblastn_fwd_1_2, tblastn_fwd_1_2_3)


tblastn_fwd_1_3 <- pair_overlaps(tblastn_fwd_1, tblastn_fwd_3)
tblastn_fwd_2_3 <- pair_overlaps(tblastn_fwd_2, tblastn_fwd_3)

tblastn_rev <- tblastn_fixed %>%
  filter(strand == "-") %>%
  as_granges()



# Operations on forward frame
tblastn_fwd_ranges <- GenomicRanges::reduce(tblastn_fwd, min.gapwidth = 10)
intersect_ranges(tblastn_fwd, tblastn_fwd_ranges)







# get filtered seq
tblastn_ranges <- GenomicRanges::reduce(tblastn_fixed)
tblastn_seq <- getSeq(genome_seq, tblastn_ranges)
names(tblastn_seq) <- paste0(seqnames(tblastn_ranges), ":", ranges(tblastn_ranges), "(", strand(tblastn_ranges), ")")