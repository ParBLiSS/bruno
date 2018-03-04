#!/bin/bash

DATA=/project/tpan7/data/human/gage_chr14/gage_human_chr14_frag_1.fastq

for t in 3 2 1
do

# run the experiments.
for exp in "_freq_minimizer" "_freq_incr_minimizer" "_freq" "" "_incr" "_freq_incr" 
do
	mpirun -np ${1} bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -C -T -L ${t} -O a4_k31${exp}_chr14_L${t}.p${1} $DATA 2>&1 | tee a4_k31${exp}_chr14_L${t}.p${1}.log
done

# compare to base case
for exp in "_freq" "_incr" 
do
	for suf in ".0.valid" ".0.valid.debug" "_branch.fasta" "_chain.components" "_chain.fasta" "_compressed_chain.debug"
	do
		diff a4_k31${exp}_chr14_L${t}.p${1}${suf} a4_k31_chr14_L${t}.p${1}${suf} > a4_k31${exp}_chr14_L${t}.p${1}${suf}.diff
	done
done
for exp in "_minimizer" "_incr"
do
#compare freq vs freq_incremental only.
for suf in ".0.valid" ".0.valid.debug" "_branch.fasta" "_branch.edges" "_chain.components" "_chain.fasta" "_chain.edges" "_compressed_chain.debug"
do
        diff a4_k31_freq_chr14_L${t}.p${1}${suf} a4_k31_freq${exp}_chr14_L${t}.p${1}${suf} > a4_k31_freq_vs${exp}_chr14_L${t}.p${1}${suf}.diff
done
done

done
