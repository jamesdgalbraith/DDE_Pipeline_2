#!/bin/bash

if [ -z "$THREADS" ]
then
    THREADS=1
fi

# make directories
mkdir -p out/split_out/ out/tblastn/

# Count contigs in genome
CONTIGS="`grep -c '>' seq/${GENOME}.fasta`"

# create database
makeblastdb -dbtype nucl -in seq/${GENOME}.fasta

# Search in parallel
echo "Initial search"
parallel --bar --jobs ${THREADS} -a data/queries.txt tblastn -query data/search_queries/{} -db seq/${GENOME}.fasta -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue frames\" -out out/split_out/{}_in_${GENOME}.out -max_target_seqs ${CONTIGS}

# combine blast output
cat out/split_out/*_in_${GENOME}.out > out/tblastn/compiled_in_${GENOME}.out
rm out/split_out/*_in_${GENOME}.out

# piece together 
echo "Getting sequences"
if [ -z "$SPECIES" ]
then
    Rscript proper_multiframe_solved_future_lapply.R -g ${GENOME} -t ${THREADS}
else
    Rscript proper_multiframe_solved_future_lapply.R -g ${GENOME} -s ${SPECIES} -t ${THREADS}
fi

# replace stop codons with ambigious
sed -i 's/\*/X/g' out/plain_tblastn_initial_fastas/${GENOME}_seq.fasta

# create dir
mkdir -p out/recip/

# blast against original query
echo "Reciprocal search"
blastp -num_threads ${THREADS} -query out/plain_tblastn_initial_fastas/${GENOME}_seq.fasta -db data/compiled_database/compiled.fasta -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen" -out out/recip/${GENOME}.out

# classify repeats
echo "Final classification"
if [ -z "$SPECIES" ]
then
    Rscript recip_search.R -g ${GENOME}
else
    Rscript recip_search.R -g ${GENOME} -s ${SPECIES}
fi

# # find orfs
# if [ -z "$SPECIES" ]
# then
#     Rscript orf_finder.R -g ${GENOME}
# else
#     Rscript orf_finder.R -g ${GENOME} -s ${SPECIES}
# fi
# echo "Complete"
