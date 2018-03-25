#!/bin/bash
#----------------------------------------------------
# Example SLURM job script to run hybrid applications
# (MPI/OpenMP or MPI/pthreads) on TACC's Stampede
# system.
#----------------------------------------------------
#SBATCH -J benchmark_bruno_qual    # Job name
##SBATCH -o bruno.o%j # Name of stdout output file(%j expands to jobId)

##SBATCH -e bruno.e%j # Name of stderr output file(%j expands to jobId)
##SBATCH -p normal    # Queue name

## OVERRIDE THESE ON COMMANDLINE.
#SBATCH -N 1               # Total number of nodes requested (16 cores/node)
#SBATCH -n 64              # Total number of mpi tasks requested

#SBATCH -t 4:00:00       # Run time (hh:mm:ss) - 1.5 hours
# The next line is required if the user has more than one project
# #SBATCH -A A-yourproject  # <-- Allocation name to charge job against




for DATASET in "S_aureus" "R_sphaeroides" 
do


ref=/project/tpan7/data/gage/${DATASET}/Data/original/genome.fasta

for run in "clean" "clean_recompact"
do

for clean in ".clean" ""
do

for L in 4 3 2 8 1 
do


contig=/scratch/tpan7/bruno/gage/${DATASET}/quality/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}.fasta
branch=/scratch/tpan7/bruno/gage/${DATASET}/quality/${DATASET}_A4_K31_L${L}_P64_freq_${run}.branch${clean}.fasta
log=/scratch/tpan7/bruno/gage/${DATASET}/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}.quast.log
outdir=/scratch/tpan7/bruno/gage/${DATASET}/quast/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}

mkdir -p $outdir


if [ ! -f ${log} ] 
  then

echo "/home/rnihalani3/graph_walking/assembly_evaluation/quast/quast.py --fast -t 64 -o ${outdir} -R ${ref} ${contig} ${branch} >> $log 2>&1"
echo "/home/rnihalani3/graph_walking/assembly_evaluation/quast/quast.py --fast -t 64 -o ${outdir} -R ${ref} ${contig} ${branch} >> $log 2>&1" > $log
/home/rnihalani3/graph_walking/assembly_evaluation/quast/quast.py --fast -t 64 -o ${outdir} -R ${ref} ${contig} ${branch} >> $log 2>&1

else
  echo "${log} exists.  skipping."
fi


done

done

done

done




for DATASET in "H_sapiens_chr14"
do


ref=/project/tpan7/data/gage/${DATASET}/Data/original/genome.fasta

for run in "clean" "clean_recompact"
do

for clean in ".clean" ""
do

for L in 1 2 3 4 8
do


contig=/scratch/tpan7/bruno/gage/${DATASET}/quality/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}.fasta
branch=/scratch/tpan7/bruno/gage/${DATASET}/quality/${DATASET}_A4_K31_L${L}_P64_freq_${run}.branch${clean}.fasta
log=/scratch/tpan7/bruno/gage/${DATASET}/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}.quast.log
outdir=/scratch/tpan7/bruno/gage/${DATASET}/quast/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}

mkdir -p $outdir

if [ ! -f ${log} ] 
  then


echo "/home/rnihalani3/graph_walking/assembly_evaluation/quast/quast.py --no-plots --no-html --no-snps --no-gc -t 64 --eukaryote -o ${outdir} -R ${ref} ${contig} ${branch} >> $log 2>&1"
echo "/home/rnihalani3/graph_walking/assembly_evaluation/quast/quast.py --no-plots --no-html --no-snps --no-gc -t 64 --eukaryote -o ${outdir} -R ${ref} ${contig} ${branch} >> $log 2>&1" > $log
/home/rnihalani3/graph_walking/assembly_evaluation/quast/quast.py --no-plots --no-html --no-snps --no-gc -t 64 --eukaryote -o ${outdir} -R ${ref} ${contig} ${branch} >> $log 2>&1

else
  echo "${log} exists.  skipping."
fi



done

done

done

done
