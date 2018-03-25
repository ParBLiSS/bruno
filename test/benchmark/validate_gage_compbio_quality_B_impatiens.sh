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
#SBATCH -n 1              # Total number of mpi tasks requested

#SBATCH -t 4:00:00       # Run time (hh:mm:ss) - 1.5 hours
# The next line is required if the user has more than one project
# #SBATCH -A A-yourproject  # <-- Allocation name to charge job against


cd /project/tpan7/data/gage/MUMmer3.23


for DATASET in "B_impatiens"
do

if [ "$DATASET" == "S_aureus" ]
then
size=2872915

elif [ "$DATASET" == "R_sphaeroides" ]
then
size=4603060

elif [ "$DATASET" == "H_sapiens_chr14" ]
then
size=88289540

elif [ "$DATASET" == "B_impatiens" ]
then
size=250000000


fi


ref=/project/tpan7/data/gage/${DATASET}/Data/original/genome.fasta

for run in "clean_incr" "clean_recompact_incr"
do

for clean in ".clean" ""
do

for L in 2 3 4 8
do


contig=/scratch/tpan7/bruno/gage/${DATASET}/quality/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}.fasta
log=/scratch/tpan7/bruno/gage/${DATASET}/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}.quality.log


echo "java GetFastaStats -o -min 200 -genomeSize $size $contig > $log 2>&1"
echo "java GetFastaStats -o -min 200 -genomeSize $size $contig > $log 2>&1" > $log
java -Xmx2048m GetFastaStats -o -min 200 -genomeSize $size $contig >> $log 2>&1


done

done

done

done
