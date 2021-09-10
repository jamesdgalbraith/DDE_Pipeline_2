#!/bin/bash

# Count contigs in genome
CONTIGS="`grep -c '>' seq/${GENOME}.fasta`"

# Search in parallel
parallel --bar --jobs 8 -a data/single_split.txt tblastn -in_pssm data/representative/psiblast_pssms/{}.smp -db seq/${GENOME}.fasta -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue frames\" -out out/pssm_out/{}_in_${GENOME}.out -max_target_seqs ${CONTIGS}

# Compile psi-tblastn out
cat out/pssm_out/*_in_${GENOME}.out > out/tblastn_pssm/compiled_in_${GENOME}.out