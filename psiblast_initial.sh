#!/bin/bash

psiblast -query data/representative/single_split/${QUERY}.fasta -db data/unaligned/compiled_database/compiled.fasta -evalue 0.001\
    -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue frames" -num_iterations 3 \
    -out_pssm data/representative/psiblast_pssms/${QUERY}.smp -out data/representative/psiblast_out/${QUERY}.out

sed -i "s/Query_1/$QUERY/" data/representative/psiblast_pssms/${QUERY}.smp
