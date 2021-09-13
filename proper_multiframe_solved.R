#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})


# parse input variables
option_list = list(
  make_option(c("-s", "--species_name"), type="character", default=NULL, 
              help="species name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out/plain_tblastn_initial_fastas/", 
              help="path to output [default= %default]", metavar="character")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.na(opt$species_name)) {
  stop("Species name is needed")
} else {
  print(opt$species_name)
}

# set species names
species_name <- opt$species_name

# read in blast output
tblastn_in <- read_tsv(paste0("out/tblastn/compiled_in_", species_name, ".out"), show_col_types = F,
                       col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue", "frames"))

# read in genome
genome_seq <- readDNAStringSet(paste0("seq/", species_name, ".fasta"))
names(genome_seq) <- gsub(" .*", "", names(genome_seq))

# filter blast output
tblastn_fixed <- tblastn_in %>%
  dplyr::filter(length >= 0.5*qlen, length <= 1.2*qlen) %>%
  tidyr::separate(frames, into = c("qframe", "sframe"), sep = "/") %>%
  dplyr::mutate(strand = ifelse(sstart < send, "+", "-"),
                start = ifelse(strand == "+", sstart, send),
                end = ifelse(strand == "+", send, sstart),
                qframe = as.integer(qframe), sframe = as.integer(sframe))

# FORWARD
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
single_frames_fwd <- tibble(joined = names(table(tblastn_fwd_combined$joined)), n = as.integer(table(tblastn_fwd_combined$joined))) %>%
  filter(n==1) %>%
  tidyr::separate(joined, into = c("seqnames", "ranges", "strand"), sep = ":") %>%
  tidyr::separate(ranges, into = c("start", "end"), sep = "-") %>%
  dplyr::mutate(start = as.integer(start), end = as.integer(end)) %>%
  plyranges::as_granges()

# Get sequence of single frames
single_frames_fwd_seq <- getSeq(genome_seq, single_frames_fwd)
single_frames_fwd_seq <- translate(single_frames_fwd_seq, if.fuzzy.codon = "X")
names(single_frames_fwd_seq) <-paste0(seqnames(single_frames_fwd), ":", ranges(single_frames_fwd), "(", strand(single_frames_fwd), ")")

# Two frames
multiple_frames_fwd <- tibble(joined = names(table(tblastn_fwd_combined$joined)), n = as.integer(table(tblastn_fwd_combined$joined))) %>%
  filter(n==2) %>%
  tidyr::separate(joined, into = c("seqnames", "ranges", "strand"), sep = ":") %>%
  tidyr::separate(ranges, into = c("start", "end"), sep = "-") %>%
  dplyr::mutate(start = as.integer(start), end = as.integer(end)) %>%
  plyranges::as_granges()

# extract seperate frames
tblastn_fwd_framed <- tblastn_fwd_combined$framed %>% mutate(frame  = tblastn_fwd_combined$frame) %>%
  sort()

# loop over, extract frames
for(j in 1:length(multiple_frames_fwd)){
  print(paste0(j, " of ", length(multiple_frames_fwd)))
  
  # first remove those wholly contained within others
  if(width(multiple_frames_fwd[j]) == width(join_overlap_intersect_directed(tblastn_fwd_framed, multiple_frames_fwd[j])[1]) ||
     width(multiple_frames_fwd[j]) == width(join_overlap_intersect_directed(tblastn_fwd_framed, multiple_frames_fwd[j])[2])){
    ship_seq <- translate(getSeq(genome_seq, multiple_frames_fwd[j]))
    names(ship_seq) <- paste0(seqnames(multiple_frames_fwd[j]), ":", ranges(multiple_frames_fwd[j]), "(", strand(multiple_frames_fwd[j]), ")")
    
  } else{
    
    # identify rear portion
    stern <- join_overlap_intersect_directed(tblastn_fwd_framed, multiple_frames_fwd[j])[2]
    
    #identify front portion and trim off (subtract) rear portion
    bow <- setdiff_ranges_directed(join_overlap_intersect_directed(tblastn_fwd_framed, multiple_frames_fwd[j])[1],
                                     join_overlap_intersect_directed(tblastn_fwd_framed, multiple_frames_fwd[j])[2])
    
    # get sequences, translate and append
    bow_seq <- translate(getSeq(genome_seq, bow))
    stern_seq <- translate(getSeq(genome_seq, stern))
    ship_seq <- AAStringSet(paste0(as.character(bow_seq), as.character(stern_seq)))
    names(ship_seq) <- paste0(seqnames(multiple_frames_fwd[j]), ":", ranges(multiple_frames_fwd[j]), "(", strand(multiple_frames_fwd[j]), ")")
    
  }
  
  
  if(j==1){multiple_frames_fwd_seq <- ship_seq}else{multiple_frames_fwd_seq <- c(multiple_frames_fwd_seq, ship_seq)}
  
}

# Compile forward sequences
fwd_seq <- c(single_frames_fwd_seq, multiple_frames_fwd_seq)

# REVERSE
# Split into separate into each frame
tblastn_rev <- tblastn_fixed %>%
  dplyr::filter(strand == "-") %>%
  plyranges::as_granges()

# Seperate into individual frames
tblastn_rev_1 <- GenomicRanges::reduce(tblastn_rev[tblastn_rev$sframe == -1]) %>% dplyr::mutate(sframe = -1)
tblastn_rev_2 <- GenomicRanges::reduce(tblastn_rev[tblastn_rev$sframe == -2]) %>% dplyr::mutate(sframe = -2)
tblastn_rev_3 <- GenomicRanges::reduce(tblastn_rev[tblastn_rev$sframe == -3]) %>% dplyr::mutate(sframe = -3)

# Reduce/merge (will contain different frames)
tblastn_rev_reduced <- GenomicRanges::reduce(tblastn_rev)

# Overlap frames and name columns
tblastn_rev_combined <- rbind(plyranges::pair_overlaps(tblastn_rev_reduced, tblastn_rev_1),
                              plyranges::pair_overlaps(tblastn_rev_reduced, tblastn_rev_2)) %>%
  rbind(plyranges::pair_overlaps(tblastn_rev_reduced, tblastn_rev_3))
colnames(tblastn_rev_combined) <- c("joined", "framed", "frame")

# Determine those in single frames and make ranges
single_frames_rev <- tibble(joined = names(table(tblastn_rev_combined$joined)), n = as.integer(table(tblastn_rev_combined$joined))) %>%
  filter(n==1) %>%
  tidyr::separate(joined, into = c("seqnames", "ranges", "strand"), sep = ":") %>%
  tidyr::separate(ranges, into = c("start", "end"), sep = "-") %>%
  dplyr::mutate(start = as.integer(start), end = as.integer(end)) %>%
  plyranges::as_granges()

# Get sequence of single frames
single_frames_rev_seq <- getSeq(genome_seq, single_frames_rev)
single_frames_rev_seq <- translate(single_frames_rev_seq, if.fuzzy.codon = "X")
names(single_frames_rev_seq) <-paste0(seqnames(single_frames_rev), ":", ranges(single_frames_rev), "(", strand(single_frames_rev), ")")

# Two frames
multiple_frames_rev <- tibble(joined = names(table(tblastn_rev_combined$joined)), n = as.integer(table(tblastn_rev_combined$joined))) %>%
  filter(n==2) %>%
  tidyr::separate(joined, into = c("seqnames", "ranges", "strand"), sep = ":") %>%
  tidyr::separate(ranges, into = c("start", "end"), sep = "-") %>%
  dplyr::mutate(start = as.integer(start), end = as.integer(end)) %>%
  plyranges::as_granges()

# extract seperate frames
tblastn_rev_framed <- tblastn_rev_combined$framed %>% mutate(frame  = tblastn_rev_combined$frame) %>%
  sort()

# loop over, extract frames
for(j in 1:length(multiple_frames_rev)){
  print(paste0(j, " of ", length(multiple_frames_rev)))
  
  # first remove those wholly contained within others
  if(width(multiple_frames_rev[j]) == width(join_overlap_intersect_directed(tblastn_rev_framed, multiple_frames_rev[j])[1]) ||
     width(multiple_frames_rev[j]) == width(join_overlap_intersect_directed(tblastn_rev_framed, multiple_frames_rev[j])[2])){
    ship_seq <- translate(getSeq(genome_seq, multiple_frames_rev[j]))
    names(ship_seq) <- paste0(seqnames(multiple_frames_rev[j]), ":", ranges(multiple_frames_rev[j]), "(", strand(multiple_frames_rev[j]), ")")
    
  } else{
    
    # identify rear portion
    stern <- join_overlap_intersect_directed(tblastn_rev_framed, multiple_frames_rev[j])[1]
    
    #identify front portion and trim off (subtract) rear portion
    bow <- setdiff_ranges_directed(join_overlap_intersect_directed(tblastn_rev_framed, multiple_frames_rev[j])[2],
                                   join_overlap_intersect_directed(tblastn_rev_framed, multiple_frames_rev[j])[1])
    
    # get sequences, translate and append
    bow_seq <- translate(getSeq(genome_seq, bow))
    stern_seq <- translate(getSeq(genome_seq, stern))
    ship_seq <- AAStringSet(paste0(as.character(bow_seq), as.character(stern_seq)))
    names(ship_seq) <- paste0(seqnames(multiple_frames_rev[j]), ":", ranges(multiple_frames_rev[j]), "(", strand(multiple_frames_rev[j]), ")")
    
  }
  
  
  if(j==1){multiple_frames_rev_seq <- ship_seq}else{multiple_frames_rev_seq <- c(multiple_frames_rev_seq, ship_seq)}
  
}

# Compile reverse sequences
rev_seq <- c(single_frames_rev_seq, multiple_frames_rev_seq)

# ALL
# Compile together
both_seq <- c(fwd_seq, rev_seq)

# Write to file
writeXStringSet(both_seq, filepath = paste0(opt$out, "", species_name, "_seq.fasta"))
