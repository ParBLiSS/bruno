#!/bin/bash

DATA=/project/tpan7/data/human/gage_chr14/gage_human_chr14_frag_1.fastq

for t in 2 
do

# run the experiments.
for exp in "_freq"
do
	mpirun -np ${1} bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -C -T -L ${t} -O a4_k31${exp}_chr14_L${t}.p${1} $DATA 2>&1 | tee a4_k31${exp}_chr14_L${t}.p${1}.log
done

# compare to base case
for exp in "_freq"
do
	for suf in ".0.valid" ".0.valid.debug" "_branch.fasta" "_chain.components" "_chain.fasta" "_compressed_chain.debug"
	do
		diff a4_k31${exp}_chr14_L${t}.p64${suf} a4_k31_chr14_L${t}.p64${suf} > a4_k31${exp}_chr14_L${t}.p64${suf}.diff
	done
done

#compare freq vs freq_incremental only.
#for suf in ".0.valid" ".0.valid.debug" "_branch.fasta" "_chain.components" "_chain.fasta" "_compressed_chain.debug"
#do
#	diff a4_k31_freq_chr14_L${t}.p64${suf} a4_k31_freq_incr_chr14_L${t}.p64${suf} > a4_k31_freq_vs_incr_chr14_L${t}.p64${suf}.diff
#done

done
