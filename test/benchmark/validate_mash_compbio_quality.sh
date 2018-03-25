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
#SBATCH -n 8              # Total number of mpi tasks requested

#SBATCH -t 4:00:00       # Run time (hh:mm:ss) - 1.5 hours
# The next line is required if the user has more than one project
# #SBATCH -A A-yourproject  # <-- Allocation name to charge job against


cd /project/tpan7/data/gage/MUMmer3.23


for DATASET in "S_aureus" "R_sphaeroides" "H_sapiens_chr14"
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

fi


ref=/project/tpan7/data/gage/${DATASET}/Data/original/genome.fasta

for run in "clean" "clean_recompact"
do

for clean in ".clean" ""
do

for L in 4 3 2 8 1
do


contig=/scratch/tpan7/bruno/gage/${DATASET}/quality/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}.fasta
branch=/scratch/tpan7/bruno/gage/${DATASET}/quality/${DATASET}_A4_K31_L${L}_P64_freq_${run}.branch${clean}.fasta
log=/scratch/tpan7/bruno/gage/${DATASET}/${DATASET}_A4_K31_L${L}_P64_freq_${run}.chain${clean}.mash.log

if [ ! -f ${log} ]
  then



echo "JACCARD WITH CONTIG AND BRANCH: $contig"
echo "JACCARD WITH CONTIG AND BRANCH" > $log
echo "/project/tpan7/mash-Linux64-v2.0/mash dist -p 8 -k 31 -g ${size} ${ref} ${contig} ${branch} >> $log 2>&1" >> $log
/project/tpan7/mash-Linux64-v2.0/mash dist -p 8 -k 31 -g ${size} ${ref} ${contig} ${branch} >> $log 2>&1

echo "JACCARD WITH CONTIG ONLY: $contig"
echo "JACCARD WITH CONTIG ONLY" >> $log
echo "/project/tpan7/mash-Linux64-v2.0/mash dist -p 8 -k 31 -g ${size} ${ref} ${contig} >> $log 2>&1" >> $log
/project/tpan7/mash-Linux64-v2.0/mash dist -p 8 -k 31 -g ${size} ${ref} ${contig} >> $log 2>&1

echo "SKETCH contig and branches: $contig"
echo "SKETCH WITH CONTIG AND BRANCH" >> $log
echo "/project/tpan7/mash-Linux64-v2.0/mash sketch -p 8 -k 31 -g ${size} ${contig} ${branch} >> $log 2>&1" >> $log
/project/tpan7/mash-Linux64-v2.0/mash sketch -p 8 -k 31 -g ${size} ${contig} ${branch} >> $log 2>&1

echo "SCREEN WITH CONTIG AND BRANCH: $contig"
echo "SCREEN WITH CONTIG AND BRANCH" >> $log
echo "/project/tpan7/mash-Linux64-v2.0/mash screen -p 8 ${contig}.msh ${ref} >> $log 2>&1" >> $log
/project/tpan7/mash-Linux64-v2.0/mash screen -p 8 ${contig}.msh ${ref} >> $log 2>&1




echo "SKETCH contig only: $contig"
echo "SKETCH WITH CONTIG only" >> $log
echo "/project/tpan7/mash-Linux64-v2.0/mash sketch -p 8 -k 31 -g ${size} ${contig} >> $log 2>&1" >> $log
/project/tpan7/mash-Linux64-v2.0/mash sketch -p 8 -k 31 -g ${size} ${contig} >> $log 2>&1

echo "SCREEN WITH CONTIG only: $contig"
echo "SCREEN WITH CONTIG only" >> $log
echo "/project/tpan7/mash-Linux64-v2.0/mash screen -p 8 ${contig}.msh ${ref} >> $log 2>&1" >> $log
/project/tpan7/mash-Linux64-v2.0/mash screen -p 8 ${contig}.msh ${ref} >> $log 2>&1


echo "SKETCH ref only: $ref"
echo "SKETCH ref only" >> $log
echo "/project/tpan7/mash-Linux64-v2.0/mash sketch -p 8 -k 31 -g ${size} ${ref} >> $log 2>&1" >> $log
/project/tpan7/mash-Linux64-v2.0/mash sketch -p 8 -k 31 -g ${size} ${ref} >> $log 2>&1

echo "SCREEN WITH ref in CONTIG AND BRANCH: $contig"
echo "SCREEN WITH ref in CONTIG AND BRANCH" >> $log
echo "/project/tpan7/mash-Linux64-v2.0/mash screen -p 8 ${ref}.msh ${contig} ${branch} >> $log 2>&1" >> $log
/project/tpan7/mash-Linux64-v2.0/mash screen -p 8 ${ref}.msh ${contig} ${branch} >> $log 2>&1

echo "SCREEN WITH ref in CONTIG : $contig"
echo "SCREEN WITH ref in CONTIG " >> $log
echo "/project/tpan7/mash-Linux64-v2.0/mash screen -p 8 ${ref}.msh ${contig} >> $log 2>&1" >> $log
/project/tpan7/mash-Linux64-v2.0/mash screen -p 8 ${ref}.msh ${contig} >> $log 2>&1



else
  echo "${log} exists.  skipping."
fi



done

done

done

done
