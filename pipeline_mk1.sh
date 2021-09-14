#!/bin/bash

# make directories
mkdir -p out/split_out/ out/tblastn/

# Count contigs in genome
CONTIGS="`grep -c '>' seq/${GENOME}.fasta`"

# create database
makeblastdb -dbtype nucl -in seq/${GENOME}.fasta

# Search in parallel
parallel --bar --jobs ${THREADS} -a data/queries.txt tblastn -query data/representative/search_queries/{} -db seq/${GENOME}.fasta -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue frames\" -out out/split_out/{}_in_${GENOME}.out -max_target_seqs ${CONTIGS} -num_threads 1

# combine blast output
cat out/split_out/*_in_${GENOME}.out > out/tblastn/compiled_in_${GENOME}.out

# piece together 
Rscript proper_multiframe_solved.R -s ${GENOME}

# replace stop codons with ambigious
sed -i 's/\*/X/g' out/plain_tblastn_initial_fastas/${GENOME}_seq.fasta

# create dir
mkdir -p out/recip/

# blast against original query
blastp -num_threads ${THREADS} -query out/plain_tblastn_initial_fastas/${GENOME}_seq.fasta -db data/representative/compiled_database/compiled.fasta -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen" -out out/recip/${GENOME}.out

# classify repeats
Rscript recip_search.R -s ${GENOME}