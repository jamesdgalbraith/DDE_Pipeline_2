#!/bin/bash

echo $SPECIES

if [ -z "$THREADS" ]
then
    THREADS=1
fi

# end if blast output is empty
if [ ! -s out/tblastn/compiled_in_${GENOME}.out ];
    then echo "No DDEs found in ${GENOME}" && exit 2
fi

# check genome present
if [ ! -s seq/${GENOME}.fasta.gz ];
    then echo "${GENOME} file is missing" && exit 2
fi

# unzip genom
gunzip < seq/${GENOME}.fasta.gz > seq/${GENOME}.fasta

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