#!/bin/bash
#PBS -l nodes=1:ppn=28
#PBS -l walltime=04:00:00
#PBS -q swarm
#PBS -n

module load gcc/4.9.4
module load hwloc
module load mvapich2/2.3b

DATAROOT=/nv/hswarm1/tpan7/scratch/bruno/data/gage
#DATA=/project/tpan7/data/human/gage_chr14/gage_human_chr14_frag_1.fastq
BINDIR=/nv/hswarm1/tpan7/data/build/bruno
OUTROOT=/nv/hswarm1/tpan7/scratch/bruno/gage

p=$PBS_NP


############### scaling - chr14 1 to 64 cores, 3 iterations, all threshold values, freq_clean_recompact, freq_clean, freq_minimizer, and baseline.
for DATASET in "H_sapiens_chr14" "B_impatiens"
do

# clear the file cache first - don't bother dropping the cache


DATA=${DATAROOT}/${DATASET}/all.fastq
OUTDIR=${OUTROOT}/${DATASET}

# warm up.
	mpirun_rsh -hostfile=$PBS_NODEFILE -np $p ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31_freq -R -B -T -L 4 -O ${OUTDIR}/${DATASET}_A4_K31_L4_P${p}_freq $DATA > ${OUTDIR}/${DATASET}_A4_K31_L4_P${p}_freq.log 2>&1


for it in 1 2 3 
do

mkdir -p ${OUTDIR}/${it}
cd ${OUTDIR}/${it}


for t in 8 4 3 2 1
do

# run the experiments.
for exp in "_freq_clean_recompact" "_freq_clean" "_freq_minimizer" "" #"_freq"
do

	echo "mpirun_rsh -hostfile=$PBS_NODEFILE -np $p ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log 2>&1"
	echo "mpirun_rsh -hostfile=$PBS_NODEFILE -np $p ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log 2>&1" > ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log
	mpirun_rsh -hostfile=$PBS_NODEFILE -np $p ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA >> ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log 2>&1
done

done

done 

done



