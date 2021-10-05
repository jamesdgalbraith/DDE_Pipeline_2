#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

# parse input variables
option_list = list(
  make_option(c("-g", "--genome_name"), type="character", default=NULL, 
              help="genome name", metavar="character"),
  make_option(c("-s", "--species_name"), type="character", default=NULL,
              help="genome name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out/orfs/", 
              help="path to output [default= %default]", metavar="character"),
  make_option(c("-f", "--flank_size"), type="integer", default=450, 
              help="Size of flanking sequence in nt (must be divisible by 3)")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
outdir <- opt$out

if (is.null(opt$genome_name)) {
  stop("Genome name is needed")
} else {
  # set genome names
  genome_name <- opt$genome_name
}

if (is.null(opt$species_name)) {
  species_name <- opt$genome_name
} else {
  species_name <- opt$species_name
}

message(species_name)

if(opt$flank_size%%3 != 0){
  stop("Size of flanking sequence in nt must be divisible by 3")
} else {
  flank <- opt$flank_size
}

if(!dir.exists(paste0(outdir, "/aa"))){
  dir.create(paste0(outdir, "/aa"), recursive = T)
}

if(!dir.exists(paste0(outdir, "/nt"))){
  dir.create(paste0(outdir, "/nt"), recursive = T)
}

# genome_name = "latCor_2.0"
# species_name = "Laticauda_colubrina"
# outdir <- "out/orfs/"

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
  library(ORFik)
})

if(!file.exists(paste0("out/plain_tblastn_initial_fastas/", genome_name, "_single_frame.bed"))){
  stop("BED contains location of transposes identified in a single frame needed")
}

# Read in genome
genome_seq <- readDNAStringSet(paste0("seq/", genome_name, ".fasta"))
names(genome_seq) <- sub(" .*", "", names(genome_seq))

# read in bed as tibble (+1 to correct start)
suppressMessages(
  single_frame_bed <- read_tsv(paste0("out/plain_tblastn_initial_fastas/", genome_name, "_single_frame.bed"),
                               col_names = c("seqnames", "start", "end", "name", "score", "strand")) %>%
    mutate(start = start + 1, name = paste0(seqnames, ":", start, "-", end, "(", strand, ")"),
           frame = case_when(start%%3 == 1 ~ 1, start%%3 == 2 ~ 2, start%%3 == 0 ~ 3)) %>%
    inner_join(tibble(seqnames = names(genome_seq), contig_width = width(genome_seq)))
  )

# read in classification of all, join to keep only singlke frames
suppressMessages(
  single_frame_classified <- tibble(name = names(readAAStringSet(paste0("out/classified_tnps/", species_name, "_compiled.fasta")))) %>%
    separate(name, into = c("name", "class"), sep = "#") %>%
    inner_join(single_frame_bed)
  )

# extend flanks by multiple of 3
single_frame_classified <- single_frame_classified %>%
  mutate(start = start - flank, end = end + flank, width = end - start +1) %>%
  filter(start > 1, end <= contig_width)

# convert classified to ranges
single_frame_classified_ranges <- as_granges(single_frame_classified)

# get sequence of hits in a single frame
single_frame_classified_seq <- getSeq(genome_seq, single_frame_classified_ranges)
names(single_frame_classified_seq) <- single_frame_classified_ranges$name

# look for orfs in sequences
classified_seq_orfs <- findORFs(seqs = single_frame_classified_seq, minimumLength = 100, longestORF = TRUE)
classified_seq_ranges <- unlist(x = classified_seq_orfs, use.names = TRUE)

# Convert IRanges to tibble and add names, ensure orfs are in reading frame 1
classified_seq_tibble <- tibble(name = names(single_frame_classified_seq)[as.integer(names(classified_seq_ranges))],
                                orf_start = start(classified_seq_ranges), orf_end = end(classified_seq_ranges),
                                orf_frame =
                                  case_when(start(classified_seq_ranges)%%3 == 1 ~ 1,
                                            start(classified_seq_ranges)%%3 == 2 ~ 2,
                                            start(classified_seq_ranges)%%3 == 0 ~ 3)) %>%
  filter(orf_frame == 1) %>%
  mutate(seqnames = sub(":.*", "", name), ranges = sub(".*:", "", name)) %>%
  separate(ranges, into = c("ranges", "strand"), sep = "\\(") %>%
  separate(ranges, into = c("start", "end"), sep = "-") %>%
  mutate(strand = sub(")", "", strand), start = as.integer(start), end = as.integer(end))

# Ensure orfs begin before and end after domain
suppressMessages(
  single_frame_classified_orfs <- single_frame_classified %>%
    dplyr::select(-seqnames, -start, -end, -strand) %>%
    inner_join(classified_seq_tibble) %>%
    filter(orf_start <= flank) %>%
    filter(orf_end > width - flank) %>%
    select(-width)
)

# single 
single_frame_classified_orf_ranges <- as_granges(single_frame_classified_orfs)
single_frame_classified_orfs_nt <- getSeq(genome_seq, single_frame_classified_orf_ranges)
names(single_frame_classified_orfs_nt) <- paste0(single_frame_classified_orf_ranges$name)

classes <- base::unique(single_frame_classified_orf_ranges$class)

compiled_aa <- readAAStringSet(paste0("out/classified_tnps/", species_name, "_compiled.fasta"))

for(i in 1:length(classes)){
  
  tmp_aa <- compiled_aa[names(compiled_aa) %in% paste0(single_frame_classified_orf_ranges$name, "#", classes[i])]
  tmp_nt <- single_frame_classified_orfs_nt[single_frame_classified_orf_ranges$class == classes[i]]
  
  if(length(tmp_nt) > 0){
    names(tmp_nt) <- paste0(species_name, "#", names(tmp_nt))
    names(tmp_aa) <- paste0(species_name, "#", sub("#.*", "", names(tmp_aa)))
    writeXStringSet(tmp_aa, paste0("out/orfs/aa/", species_name, "_", classes[i], ".fasta"))
    writeXStringSet(tmp_nt, paste0("out/orfs/nt/", species_name, "_", classes[i], ".fasta"))
  }
  
}
