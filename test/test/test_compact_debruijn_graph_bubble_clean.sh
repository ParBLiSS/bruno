#!/bin/bash


DATA=~/src/bruno/test/data/debruijn.bubble2.fastq
for t in 4 3 2 1
do

# run the experiments.
for exp in "freq" "freq_clean" "freq_clean_recompact"
do
	echo "mpirun -np ${1} bin/compact_debruijn_graph_fastq_A4_K31_${exp} -R -B -C -T -L ${t} -O a4_k31_chr14_bub2_L${t}.p${1}.${exp} $DATA > a4_k31_chr14_bub2_L${t}.p${1}.${exp}.log 2>&1"
	mpirun -np ${1} bin/compact_debruijn_graph_fastq_A4_K31_${exp} -R -B -C -T -L ${t} -O a4_k31_chr14_bub2_L${t}.p${1}.${exp} $DATA > a4_k31_chr14_bub2_L${t}.p${1}.${exp}.log 2>&1
done

# compare to base case
for exp in "freq_clean" "freq_clean_recompact"
do
	for suf in "0.valid" "branch.fasta" "branch.edges" "chain.components" "chain.fasta" "chain.edges" "debug_compressed_chain"
	do
		diff a4_k31_chr14_bub2_L${t}.p${1}.${exp}.${suf} a4_k31_chr14_bub2_L${t}.p${1}.freq.${suf} > a4_k31_chr14_bub2_L${t}.p${1}.freq_vs_${exp}.${suf}.diff
	done
done

done
