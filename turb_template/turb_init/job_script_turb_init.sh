#!/bin/sh

#SBATCH -J !!JOB_NAME!!
#SBATCH -o !!JOB_NAME!!."%j".out
#SBATCH -e !!JOB_NAME!!."%j".err
#SBATCH --mail-user !!USER_EMAIL!!
#SBATCH --partition=!!QUEUE!!
#SBATCH --mail-type=ALL
#SBATCH --nodes=!!NODES!!
#SBATCH --ntasks-per-node=!!NTASKS_PER_NODE!!
#SBATCH --time=!!TIME_LIMIT!!

set -e
SECONDS=0

module purge
module load intel/19.1.2
module load impi/2019.8
module load fftw-mpi/3.3.8
module load hdf5-mpi/1.8.21
module list

cd !!WORK_DIR!! 

srun ./athena_turb -i !!INPUT_FILE!! -t !!TIME_LIMIT_RST!!

echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo "Boom!"
