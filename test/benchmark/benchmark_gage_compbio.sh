#!/bin/bash
#----------------------------------------------------
# Example SLURM job script to run hybrid applications 
# (MPI/OpenMP or MPI/pthreads) on TACC's Stampede 
# system.
#----------------------------------------------------
#SBATCH -J benchmark_bruno_speed    # Job name
##SBATCH -o bruno.o%j # Name of stdout output file(%j expands to jobId)
 
##SBATCH -e bruno.e%j # Name of stderr output file(%j expands to jobId)
##SBATCH -p normal    # Queue name

## OVERRIDE THESE ON COMMANDLINE.
#SBATCH -N 1               # Total number of nodes requested (16 cores/node)
#SBATCH -n 64              # Total number of mpi tasks requested

#SBATCH -t 72:00:00       # Run time (hh:mm:ss) - 1.5 hours
# The next line is required if the user has more than one project
# #SBATCH -A A-yourproject  # <-- Allocation name to charge job against


DATAROOT=/project/tpan7/data/gage
BINDIR=/nethome/tpan7/build/bruno
OUTROOT=/scratch/tpan7/bruno/gage

############### quality - all datasets, 64 cores, 1 iteration, all threshold values, freq_clean or freq_clean_recompact only.
for DATASET in "S_aureus" "R_sphaeroides" "H_sapiens_chr14"
do

# clear the file cache first
drop_caches


DATA=${DATAROOT}/${DATASET}/Data/original/all.fastq
OUTDIR=${OUTROOT}/${DATASET}

for p in 64
do

ppn=16

for it in "quality" 
do

mkdir -p ${OUTDIR}/${it}
cd ${OUTDIR}/${it}


for t in 8 4 3 2 1
do

# run the experiments.
for exp in "_freq_clean_recompact" "_freq_clean" 
do

	echo "mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log 2>&1"
	echo "mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log 2>&1" > ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log
	mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA >> ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log 2>&1
done

done

done 

done

done

############### scaling - chr14 1 to 64 cores, 3 iterations, all threshold values, freq_clean_recompact, freq_clean, freq_minimizer, and baseline.
for DATASET in "H_sapiens_chr14"
do

# data should already be in file cache from prev iterations


DATA=${DATAROOT}/${DATASET}/Data/original/all.fastq
OUTDIR=${OUTROOT}/${DATASET}

for p in 64 32 16 8 4
do

ppn=$((p / 4))


for it in 1 2 3
do

mkdir -p ${OUTDIR}/${it}
cd ${OUTDIR}/${it}


for t in 8 4 3 2 1
do

# run the experiments.
for exp in "_freq_clean_recompact" "_freq_clean" "_freq_minimizer" "" #"_freq"
do
	
	echo "mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${DATASET}_A4_K31_L${t}_P${p}${exp}.${it}.log 2>&1"
	echo "mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${DATASET}_A4_K31_L${t}_P${p}${exp}.${it}.log 2>&1" > ${OUTDIR}/${DATASET}_A4_K31_L${t}_P${p}${exp}.${it}.log
	mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA >> ${OUTDIR}/${DATASET}_A4_K31_L${t}_P${p}${exp}.${it}.log 2>&1
	
	rm ${OUTDIR}/${it}/*
done

done

done 

done

#  ppn is 1 for 2 and 1 processes.
for p in 2 1
do

ppn=1


for it in 1 2 3
do

mkdir -p ${OUTDIR}/${it}
cd ${OUTDIR}/${it}


for t in 8 4 3 2 1
do

# run the experiments.
for exp in "_freq_clean_recompact" "_freq_clean" "_freq_minimizer" "" #"_freq"
do
		
	echo "mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${DATASET}_A4_K31_L${t}_P${p}${exp}.${it}.log 2>&1"
	echo "mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${DATASET}_A4_K31_L${t}_P${p}${exp}.${it}.log 2>&1" > ${OUTDIR}/${DATASET}_A4_K31_L${t}_P${p}${exp}.${it}.log
	mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA >> ${OUTDIR}/${DATASET}_A4_K31_L${t}_P${p}${exp}.${it}.log 2>&1
	
	rm ${OUTDIR}/${it}/*
	
done

done

done 

done


done


############### quality - all datasets, 64 cores, 1 iteration, all threshold values, freq_clean or freq_clean_recompact only.
for DATASET in "B_impatiens"
do

# clear the file cache first
drop_caches


DATA=${DATAROOT}/${DATASET}/Data/original/all.fastq
OUTDIR=${OUTROOT}/${DATASET}

for p in 64
do

ppn=16

for it in "quality" 
do

mkdir -p ${OUTDIR}/${it}
cd ${OUTDIR}/${it}


for t in 8 4 3 2 1
do

# run the experiments.
for exp in "_freq_clean_recompact" "_freq_clean" 
do
	echo "mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log 2>&1"
	echo "mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA > ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log 2>&1" > ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log
	mpirun -np ${p} --map-by ppr:${ppn}:socket --rank-by core --bind-to core ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp} $DATA >> ${OUTDIR}/${it}/${DATASET}_A4_K31_L${t}_P${p}${exp}.log 2>&1

done

done

done 

done

done



