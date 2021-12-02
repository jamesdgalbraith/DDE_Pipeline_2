#!/bin/bash

echo $SPECIES

if [ -z "$THREADS" ]
then
    THREADS=1
fi

# check genome present
if [ ! -s seq/${GENOME}.fasta.gz ];
    then echo "${GENOME} file is missing" && exit 2
fi

# unzip genome
gunzip < seq/${GENOME}.fasta.gz > seq/${GENOME}.fasta

# make directories
mkdir -p out/split_out/ out/tblastn/

# Count contigs in genome
CONTIGS="`grep -c '>' seq/${GENOME}.fasta`"

# create database
makeblastdb -dbtype nucl -in seq/${GENOME}.fasta -logfile makeblast.log

# Search in parallel
echo "Initial search"
parallel --bar --jobs ${THREADS} -a data/queries.txt tblastn -query data/search_queries/{} -db seq/${GENOME}.fasta -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue frames\" -out out/split_out/{}_in_${GENOME}.out -max_target_seqs ${CONTIGS}

# combine blast output
cat out/split_out/*_in_${GENOME}.out > out/tblastn/compiled_in_${GENOME}.out
rm out/split_out/*_in_${GENOME}.out

# end if blast output is empty
if [ ! -s out/tblastn/compiled_in_${GENOME}.out ];
    then echo "No DDEs found in ${GENOME}" && exit 2
fi

# piece together 
echo "Getting sequences"
if [ -z "$SPECIES" ]
then
    Rscript proper_multiframe_solved_future_lapply.R -g ${GENOME} -t ${THREADS}
else
    Rscript proper_multiframe_solved_future_lapply.R -g ${GENOME} -s ${SPECIES} -t ${THREADS}
fi

# end if blast output is empty
if [ ! -s out/plain_tblastn_initial_fastas/${GENOME}_seq.fasta ];
    then echo "No DDEs found in ${GENOME} during stitching" && exit 2
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

# remove working seq
rm seq/${GENOME}*.nhr seq/${GENOME}*.nsq seq/${GENOME}*.nin seq/${GENOME}.fasta
