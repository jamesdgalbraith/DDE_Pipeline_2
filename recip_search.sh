#!/bin/bash

# replace stop codons with ambigious
sed -i 's/\*/X/g' out/plain_tblastn_initial_fastas/${GENOME}_seq.fasta

# create dir
mkdir -p out/recip/

# blast against original query
blastp -num_threads ${THREADS} -query out/plain_tblastn_initial_fastas/${GENOME}_seq.fasta -db data/compiled_database/compiled.fasta -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue" -out out/recip/${GENOME}.out