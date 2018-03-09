#!/bin/bash
#----------------------------------------------------
# Example SLURM job script to run hybrid applications 
# (MPI/OpenMP or MPI/pthreads) on TACC's Stampede 
# system.
#----------------------------------------------------
#SBATCH -J benchmark_bruno     # Job name
##SBATCH -o bruno.o%j # Name of stdout output file(%j expands to jobId)
 
##SBATCH -e bruno.e%j # Name of stderr output file(%j expands to jobId)
##SBATCH -p normal    # Queue name

## OVERRIDE THESE ON COMMANDLINE.
#SBATCH -N 1               # Total number of nodes requested (16 cores/node)
#SBATCH -n 64              # Total number of mpi tasks requested

#SBATCH -t 06:00:00       # Run time (hh:mm:ss) - 1.5 hours
# The next line is required if the user has more than one project
# #SBATCH -A A-yourproject  # <-- Allocation name to charge job against


DATA=/project/tpan7/data/human/gage_chr14/gage_human_chr14.fastq
BINDIR=/nethome/tpan7/build/bruno
OUTDIR=/scratch/tpan7/bruno

p=64

for it in 1 2 3
do

mkdir -p ${OUTDIR}/${it}

cd ${OUTDIR}/${it}

for t in 4 3 2 1
do

# run the experiments.
for exp in "_freq_clean_recompact" "_freq_clean" "_freq_minimizer" "_freq" "" #"_incr" "_freq_incr" "_freq_incr_minimizer" 
do

	echo "mpirun -np ${p} ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/a4_k31_chr14_L${t}.p${p}${exp} $DATA > ${OUTDIR}/${it}/a4_k31_chr14_L${t}.p${p}${exp}.log 2>&1"
	mpirun -np ${p} ${BINDIR}/bin/compact_debruijn_graph_fastq_A4_K31${exp} -R -B -T -L ${t} -O ${OUTDIR}/${it}/a4_k31_chr14_L${t}.p${p}${exp} $DATA > ${OUTDIR}/${it}/a4_k31_chr14_L${t}.p${p}${exp}.log 2>&1
done


if [ 0 -eq 1 ]
then

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

fi

done 

done
