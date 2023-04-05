#!/bin/bash -l
#SBATCH --job-name="alpha"
#SBATCH --time=16:00:00
#SBATCH --nodes=25
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --hint=nomultithread

#SBATCH -e alpha_error_file.e
#SBATCH -o alpha_output_file.o

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

for i in $(seq 54 2 76)
do
   srun python3 ./alpha_mpi.py $i
done
