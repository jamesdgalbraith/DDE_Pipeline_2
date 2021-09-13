#!/bin/bash

# replace stop codons with ambigious
sed -i 's/\*/X/g' out/plain_tblastn_initial_fastas/${GENOME}_seq.fasta

# blast against original query
blastp -num_threads 8 -query out/plain_tblastn_initial_fastas/${GENOME}_seq.fasta -db data/representative/compiled_database/compiled.fasta -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen" -out out/recip/${GENOME}.out