#!/bin/bash


DATA=~/src/bruno/test/data/debruijn.palindrome_chr14.fastq
for t in 2 1
do

# run the experiments.
for exp in "_freq" "" "_incr" "_freq_incr"
do
	mpirun -np ${1} bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -C -T -L ${t} -O a4_k31${exp}_chr14_pal_L${t}.p${1} $DATA 2>&1 | tee a4_k31${exp}_chr14_pal_L${t}.p${1}.log
done

# compare to base case
for exp in "_freq" "_incr" "_freq_incr"
do
	for suf in ".0.valid" ".0.valid.debug" "_branch.fasta" "_branch.edges" "_chain.components" "_chain.fasta" "_chain.edges" "_compressed_chain.debug"
	do
		diff a4_k31${exp}_chr14_pal_L${t}.p${1}${suf} a4_k31_chr14_pal_L${t}.p${1}${suf} > a4_k31${exp}_chr14_pal_L${t}.p${1}${suf}.diff
	done
done

#compare freq vs freq_incremental only.
for suf in ".0.valid" ".0.valid.debug" "_branch.fasta" "_branch.edges" "_chain.components" "_chain.fasta" "_chain.edges" "_compressed_chain.debug"
do
	diff a4_k31_freq_chr14_pal_L${t}.p${1}${suf} a4_k31_freq_incr_chr14_pal_L${t}.p${1}${suf} > a4_k31_freq_vs_incr_chr14_pal_L${t}.p${1}${suf}.diff
done

done
