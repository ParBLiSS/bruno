#!/bin/bash

DATA=/project/tpan7/data/human/gage_chr14/gage_human_chr14_frag_1.fastq
BINDIR=/nethome/tpan7/build/bruno
OUTDIR=/scratch/tpan7/bruno

cd ${OUTDIR}


for t in 4 #3 2 1
do

# run the experiments.
for exp in "_freq_clean_recompact" #"_freq_clean" "_freq_minimizer" "_freq" "" #"_incr" "_freq_incr" "_freq_incr_minimizer" 
do

	echo "mpirun -np ${1} ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -C -T -L ${t} -O ${OUTDIR}/a4_k31_chr14_L${t}.p${1}${exp} $DATA > ${OUTDIR}/a4_k31_chr14_L${t}.p${1}${exp}.log 2>&1"
	mpirun -np ${1} ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -C -T -L ${t} -O ${OUTDIR}/a4_k31_chr14_L${t}.p${1}${exp} $DATA > ${OUTDIR}/a4_k31_chr14_L${t}.p${1}${exp}.log 2>&1
done

# compare to base case
for exp in "_freq" #"_incr" 
do
	for suf in ".0.valid" "_branch.fasta" "_chain.components" "_chain.fasta" #"_compressed_chain.debug"
	do
		diff a4_k31_chr14_L${t}.p${1}${exp}${suf} a4_k31_chr14_L${t}.p${1}${suf} > a4_k31_chr14_L${t}.p${1}${exp}${suf}.diff
	done
done
for exp in "_minimizer" "_clean" #"_incr"
do
#compare freq vs freq_incremental only.
for suf in ".0.valid" ".0.valid.debug" "_branch.fasta" "_branch.edges" "_chain.components" "_chain.fasta" "_chain.edges" "_compressed_chain.debug"
do
        diff a4_k31_chr14_L${t}.p${1}_freq${suf} a4_k31_chr14_L${t}.p${1}_freq${exp}${suf} > a4_k31_chr14_L${t}.p${1}_freq_vs${exp}${suf}.diff
done
done

for exp in "_recompact" 
do
#compare freq vs freq_incremental only.
for suf in ".branch.fasta" ".branch.edges" ".chain.components" ".chain.fasta" ".chain.edges" #"_compressed_chain.debug"
do
        diff a4_k31_chr14_L${t}.p${1}_freq_clean${suf} a4_k31_chr14_L${t}.p${1}_freq_clean${exp}${suf} > a4_k31_chr14_L${t}.p${1}_freq_clean_vs${exp}${suf}.diff
done
done


done
