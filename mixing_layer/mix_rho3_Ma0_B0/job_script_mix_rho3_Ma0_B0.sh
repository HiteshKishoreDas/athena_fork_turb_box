#!/bin/sh
#
#SBATCH -J m300
#SBATCH -o m300."%j".out
#SBATCH -e m300."%j".err
#SBATCH --mail-user hitesh@mpa-garching.mpg.de
#SBATCH --partition=p.24h
#SBATCH --mail-type=ALL
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=23:55:00

set -e
SECONDS=0

module purge
module load intel/19.1.2
module load impi/2019.8
module load hdf5-mpi/1.8.21
module list

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_HOME/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FFTW_HOME/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$I_MPI_ROOT/intel64/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$I_MPI_ROOT/intel64/lib/release/

cd /ptmp/mpa/hitesh/athena_fork_turb_box/mixing_layer/mix_rho3_Ma0_B0/

srun ./athena_mix -i athinput_mix_rho3_Ma0_B0 -t 23:52:00
# srun ./athena -r Turb.final.rst -t 00:25:00

echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo "Boom!"
