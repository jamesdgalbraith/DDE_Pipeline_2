#!/bin/bash

# Count contigs in genome
CONTIGS="`grep -c '>' seq/${GENOME}.fasta`"

# Search in parallel
parallel --bar --jobs 4 -a data/queries.txt tblastn -query data/representative/search_queries/{}.fasta -db seq/${GENOME}.fasta -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue\" -out out/split_out/{.}_in_${GENOME}.out -max_target_seqs ${CONTIGS} -num_threads 2

cat out/split_out/*_in_${GENOME}.out > out/split_out/compiled_in_${GENOME}.out