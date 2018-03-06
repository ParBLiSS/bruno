#!/bin/bash


DATA=~/src/bruno/test/data/debruijn.bubble2.fastq
for t in 4 3 2 1
do

# run the experiments.
for exp in "_freq" "_freq_clean" "_freq_clean_recompact"
do
	mpirun -np ${1} bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -C -T -L ${t} -O a4_k31${exp}_chr14_bub2_L${t}.p${1} $DATA 2>&1 | tee a4_k31${exp}_chr14_bub2_L${t}.p${1}.log
done

# compare to base case
for exp in "_freq_clean" "_freq_clean_recompact"
do
	for suf in ".0.valid" ".0.valid.debug" "_branch.fasta" "_branch.edges" "_chain.components" "_chain.fasta" "_chain.edges" "_compressed_chain.debug"
	do
		diff a4_k31_freq${exp}_chr14_bub2_L${t}.p${1}${suf} a4_k31_freq_chr14_bub2_L${t}.p${1}${suf} > a4_k31_freq${exp}_chr14_bub2_L${t}.p${1}${suf}.diff
	done
done

done
