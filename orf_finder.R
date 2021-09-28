suppressPackageStartupMessages({
  library(optparse)
})

# genome_name = "latCor_2.0"
# species_name = "Laticauda_colubrina"
# outdir <- "out/plain_tblastn_initial_fastas/"

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
  # set genome names
  genome_name <- opt$genome_name
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
  library(ORFik)
})

if(!file.exists(paste0(outdir, "/", genome_name, "_single_frame.bed"))){
  stop("BED contains location of transposes identified in a single frame needed")
}

# Read in genome
genome_seq <- readDNAStringSet(paste0("seq/", genome_name, ".fasta"))
names(genome_seq) <- sub(" .*", "", names(genome_seq))

# read in bed as tibble (+1 to correct start)
single_frame_bed <- read_tsv(paste0("out/plain_tblastn_initial_fastas/", species_name, "_single_frame.bed"),
                             col_names = c("seqnames", "start", "end", "name", "score", "strand")) %>%
  mutate(start = start + 1, name = paste0(seqnames, ":", start, "-", end, "(", strand, ")"),
         frame = case_when(start%%3 == 1 ~ 1, start%%3 == 2 ~ 2, start%%3 == 0 ~ 3)) %>%
  inner_join(tibble(seqnames = names(genome_seq), contig_width = width(genome_seq)))

# read in classification of all, join to keep only singlke frames
single_frame_classified <- tibble(name = names(readAAStringSet(paste0("out/classified_tnps/", species_name, "_compiled.fasta")))) %>%
  separate(name, into = c("name", "class"), sep = "#") %>%
  inner_join(single_frame_bed) %>%
  mutate(start = start - 450, end = end + 450, width = end - start +1) %>%
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
single_frame_classified_orfs <- single_frame_classified %>%
  inner_join(classified_seq_tibble) %>%
  filter(orf_start <= 450) %>%
  filter(orf_end > width - 450) %>%
  arrange(orf_end - orf_start)

# single 
single_frame_classified_orf_ranges <- as_granges(classified_seq_tibble)
single_frame_classified_orf_ranges <- single_frame_classified_orf_ranges[single_frame_classified_orf_ranges$name
                                                                         %in% single_frame_classified_orfs$name]

single_frame_classified_orfs_nt <- getSeq(genome_seq, single_frame_classified_orf_ranges)
names(single_frame_classified_orfs_nt) <- paste0(single_frame_classified_orf_ranges$name)

classes <- base::unique(single_frame_classified_ranges$class)

for(i in 1:length(classes)){

  tmp <- readAAStringSet(paste0("out/classified_tnps/", classes[i], "/Laticauda_colubrina_", classes[i], ".fasta"))
  tmp_nt <- single_frame_classified_orfs_nt[paste0(single_frame_classified_orf_ranges$name, "#", classes[i]) %in% names(tmp)]
  if(length(tmp_nt) > 0){
    names(tmp_nt) <- paste0(names(tmp_nt), "#", classes[i])
    writeXStringSet(tmp_nt, paste0("out/classified_tnps/", classes[i], "/Laticauda_colubrina_", classes[i], "_nt.fasta"))
  }
  
}

clstr <- read.csv("out/classified_tnps/Tc1marPlm/Laticauda_colubrina_Tc1marPlm_nt_cd-hit-est.fasta.clstr", sep = "\t",
         row.names = NULL, header = FALSE, stringsAsFactors = FALSE) %>% as_tibble()

clstr2 <- clstr
n = nrow(clstr)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(clstr2[row,1]) == TRUE) {
    clstr2[row,1] <- x}
  else {NULL}
  x <- clstr2[row,1]
}

clstr2$V1 <- sub(">", "", clstr2$V1)
