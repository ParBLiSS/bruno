#!/bin/bash
#PBS -l nodes=1:ppn=28
#PBS -l walltime=04:00:00
#PBS -q swarm
#PBS -n

module load gcc/4.9.4
module load hwloc
module load mvapich2/2.3b


DATA=/nv/hswarm1/tpan7/data/data/gage-chr14/gage_human_chr14_frag_1.fastq
#DATA=/project/tpan7/data/human/gage_chr14/gage_human_chr14_frag_1.fastq
BINDIR=/nv/hswarm1/tpan7/data/build/bruno
OUTDIR=/nv/hswarm1/tpan7/scratch/bruno

p=$PBS_NP

cd $OUTDIR

for t in 4 3 2 1
do

# run the experiments.
for exp in "_freq_clean" "_freq_clean_recompact" "_freq_minimizer" "_freq" "" 
do
	echo "mpirun -hostfile=$PBS_NODEFILE -np $p ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -C -T -L ${t} -O ${OUTDIR}/a4_k31${exp}_chr14_L${t}.p${p} $DATA > ${OUTDIR}/a4_k31${exp}_chr14_L${t}.p${p}.log 2>&1"
	mpirun_rsh -hostfile=$PBS_NODEFILE -np $p ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -C -T -L ${t} -O ${OUTDIR}/a4_k31${exp}_chr14_L${t}.p${p} $DATA > ${OUTDIR}/a4_k31${exp}_chr14_L${t}.p${p}.log 2>&1
done

# compare to base case
for exp in "_freq" 
do
	for suf in "_branch.fasta" "_chain.fasta" 
	do
		diff a4_k31${exp}_chr14_L${t}.p${p}${suf} a4_k31_chr14_L${t}.p${p}${suf} > a4_k31${exp}_chr14_L${t}.p${p}${suf}.diff
	done
done

for exp in "_minimizer"
do
	for suf in "_branch.fasta" "_chain.fasta" 
	do
		diff a4_k31_freq${exp}_chr14_L${t}.p${p}${suf} a4_k31_freq_chr14_L${t}.p${p}${suf} > a4_k31_freq_vs${exp}_chr14_L${t}.p${p}${suf}.diff
	done
done

for exp in "_recompact"
do
	for suf in ".branch.clean.fasta" ".chain.clean.fasta" 
	do
		diff a4_k31_freq_clean${exp}_chr14_L${t}.p${p}${suf} a4_k31_freq_clean_chr14_L${t}.p${p}${suf} > a4_k31_freq_clean_vs${exp}_chr14_L${t}.p${p}${suf}.diff
	done
done



done
