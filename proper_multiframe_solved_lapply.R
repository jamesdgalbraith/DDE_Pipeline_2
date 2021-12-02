#!/usr/bin/env Rscript

ptm <- proc.time()

suppressPackageStartupMessages({
  library(optparse)
})


# parse input variables
option_list = list(
  make_option(c("-g", "--genome_name"), type="character", default=NULL,
              help="genome name", metavar="character"),
  make_option(c("-s", "--species_name"), type="character", default=NULL,
              help="genome name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out/plain_tblastn_initial_fastas/",
              help="path to output [default= %default]", metavar="character")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
outdir <- opt$out

if (is.null(opt$genome_name)) {
  stop("Genome name is needed")
} else {
  genome_name <- opt$genome_name
  
}

if (!dir.exists(outdir)){
  dir.create(outdir)
  
}

if (is.null(opt$species_name)) {
  species_name <- opt$genome_name
} else {
  species_name <- opt$species_name
}

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

# read in blast output
tblastn_fixed <- read_tsv(paste0("out/tblastn/compiled_in_", genome_name, ".out"), show_col_types = F,
                          col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend", "qlen",
                                        "sstart", "send", "slen", "evalue", "frames")) %>%
  dplyr::filter(length >= 0.5*qlen, length <= 1.2*qlen, pident >=50) %>%
  tidyr::separate(frames, into = c("qframe", "sframe"), sep = "/") %>%
  dplyr::mutate(strand = ifelse(sstart < send, "+", "-"),
                start = ifelse(strand == "+", sstart, send),
                end = ifelse(strand == "+", send, sstart),
                qframe = as.integer(qframe), sframe = as.integer(sframe))

# Split into separate into each frame
tblastn_fwd <- tblastn_fixed %>%
  dplyr::filter(strand == "+") %>%
  plyranges::as_granges()
tblastn_rev <- tblastn_fixed %>%
  dplyr::filter(strand == "-") %>%
  plyranges::as_granges()

tblastn_fixed <- NULL
invisible(gc())

# read in genome
genome_seq <- readDNAStringSet(paste0("seq/", genome_name, ".fasta"))
names(genome_seq) <- gsub(" .*", "", names(genome_seq))

# FORWARD
if(length(tblastn_fwd) > 0){
# Seperate into individual frames
tblastn_fwd_1 <- GenomicRanges::reduce(tblastn_fwd[tblastn_fwd$sframe == 1])
if(length(tblastn_fwd_1) > 0){tblastn_fwd_1$sframe <- 1}
tblastn_fwd_2 <- GenomicRanges::reduce(tblastn_fwd[tblastn_fwd$sframe == 2])
if(length(tblastn_fwd_2) > 0){tblastn_fwd_2$sframe <- 2}
tblastn_fwd_3 <- GenomicRanges::reduce(tblastn_fwd[tblastn_fwd$sframe == 3])
if(length(tblastn_fwd_3) > 0){tblastn_fwd_3$sframe <- 3}

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
if(length(multiple_frames_fwd) >0){
multiple_frames_fwd_seq <- lapply(seq_along(multiple_frames_fwd), function(j){
  # first remove those wholly contained within others
  x <- multiple_frames_fwd[j]
  
  if(width(x) == width(join_overlap_intersect_directed(tblastn_fwd_framed, x)[1]) ||
     width(x) == width(join_overlap_intersect_directed(tblastn_fwd_framed, x)[2])){
    tmp <- translate(getSeq(genome_seq, x), if.fuzzy.codon = "X")
    names(tmp) <- paste0(seqnames(x), ":", ranges(x), "(", strand(x), ")")
    
  } else{
    
    # identify rear portion
    stern <- join_overlap_intersect_directed(tblastn_fwd_framed, x)[2]
    
    #identify front portion and trim off (subtract) rear portion
    bow <- setdiff_ranges_directed(join_overlap_intersect_directed(tblastn_fwd_framed, x)[1],
                                   join_overlap_intersect_directed(tblastn_fwd_framed, x)[2])
    
    # get sequences, translate and append
    bow_seq <- suppressWarnings(translate(getSeq(genome_seq, bow), if.fuzzy.codon = "X"))
    stern_seq <- translate(getSeq(genome_seq, stern), if.fuzzy.codon = "X")
    tmp <- AAStringSet(paste0(as.character(bow_seq), as.character(stern_seq)))
    names(tmp) <- paste0(seqnames(x), ":", ranges(x), "(", strand(multiple_frames_fwd[j]), ")")
    
  }
  
  tmp
  
})

multiple_frames_fwd_seq <- do.call(c, multiple_frames_fwd_seq)

} else {
  multiple_frames_fwd_seq <- NULL
}
} else {
  
  # if neither single or multiple frames create NULL 
  single_frames_fwd_seq <- NULL
  multiple_frames_fwd_seq <- NULL
  
}

# REVERSE
if(length(tblastn_rev) > 0){
# Seperate into individual frames
tblastn_rev_1 <- GenomicRanges::reduce(tblastn_rev[tblastn_rev$sframe == -1])
if(length(tblastn_rev_1) > 0){tblastn_rev_1$sframe <- -1}
tblastn_rev_2 <- GenomicRanges::reduce(tblastn_rev[tblastn_rev$sframe == -2])
if(length(tblastn_rev_2) > 0){tblastn_rev_2$sframe <- -4}
tblastn_rev_3 <- GenomicRanges::reduce(tblastn_rev[tblastn_rev$sframe == -3])
if(length(tblastn_rev_3) > 0){tblastn_rev_3$sframe <- -3}

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
if(length(multiple_frames_rev) > 0){
multiple_frames_rev_seq <- lapply(seq_along(multiple_frames_rev), function(j){
  # first remove those wholly contained within others
  x <- multiple_frames_rev[j]

  # first remove those wholly contained within others
  if(width(x) == width(join_overlap_intersect_directed(tblastn_rev_framed, x)[1]) ||
     width(x) == width(join_overlap_intersect_directed(tblastn_rev_framed, x)[2])){
    tmp <- translate(getSeq(genome_seq, x), if.fuzzy.codon = "X")
    names(tmp) <- paste0(seqnames(x), ":", ranges(x), "(", strand(x), ")")
    
  } else{
    
    # identify rear portion
    stern <- join_overlap_intersect_directed(tblastn_rev_framed, x)[1]
    
    #identify front portion and trim off (subtract) rear portion
    bow <- setdiff_ranges_directed(join_overlap_intersect_directed(tblastn_rev_framed, x)[2],
                                   join_overlap_intersect_directed(tblastn_rev_framed, x)[1])
    
    # get sequences, translate and append
    bow_seq <- suppressWarnings(translate(getSeq(genome_seq, bow), if.fuzzy.codon = "X"))
    stern_seq <- translate(getSeq(genome_seq, stern), if.fuzzy.codon = "X")
    tmp <- AAStringSet(paste0(as.character(bow_seq), as.character(stern_seq)))
    names(tmp) <- paste0(seqnames(x), ":", ranges(x), "(", strand(x), ")")
    
  }
  
  
  tmp
  
})

multiple_frames_rev_seq <- do.call(c, multiple_frames_rev_seq)

} else {
  
  multiple_frames_rev_seq <- NULL
  
}
} else {
  
  # if neither single or multiple frames create NULL 
  single_frames_rev_seq <- NULL
  multiple_frames_rev_seq <- NULL
  
}

# Compile together
both_seq <- c(single_frames_fwd_seq, multiple_frames_fwd_seq, single_frames_rev_seq, multiple_frames_rev_seq)


# Kill if none were found
if(length(both_seq) == 0){
  
  stop("No DDEs found during stitching")
  
}

# Write to file
writeXStringSet(both_seq, filepath = paste0(outdir, "", genome_name, "_seq.fasta"))

message(paste0(as.double(proc.time() - ptm)[3], " seconds"))

# write single frame seq to file if present
if(suppressWarnings(length(c(single_frames_fwd, single_frames_rev))) > 0){
  
  write_bed(suppressWarnings(c(single_frames_fwd, single_frames_rev)) %>% select(-n),
            file = paste0(outdir, "", genome_name, "_single_frame.bed"))
  
  writeXStringSet(c(single_frames_fwd_seq, single_frames_rev_seq),
                  filepath = paste0(outdir, "", genome_name, "_single_frame.fasta"))
  
}
