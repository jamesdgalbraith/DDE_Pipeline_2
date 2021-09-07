#!/bin/bash

# Count contigs in genome
CONTIGS="`grep -c '>' seq/${GENOME}.fasta`"

# Search in parallel
parallel --bar --jobs 4 -a families.txt tblastn -query data/representative/{}.fasta -db seq/${GENOME}.fasta -evalue 1e-5 -outfmt \"6 qseqid sseqid length qstart qend qlen sstart send slen evalue\" -out out/{.}_in_${GENOME}.out -max_target_seqs ${CONTIGS} -num_threads 2