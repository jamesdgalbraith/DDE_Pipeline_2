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
#rm out/split_out/*_in_${GENOME}.out

# remove working seq
rm seq/${GENOME}*.nhr seq/${GENOME}*.nsq seq/${GENOME}*.nin seq/${GENOME}.fasta