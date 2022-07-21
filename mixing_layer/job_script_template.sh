#!/bin/sh
#
#SBATCH -J mix
#SBATCH -e mix."%j".err
#SBATCH --mail-user hitesh@mpa-garching.mpg.de
#SBATCH --partition=p.24h
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00

set -e
SECONDS=0

module purge
module load intel/19.1.2
module load impi/2019.8
module load fftw-mpi/3.3.8
module load hdf5-mpi/1.8.21
module list

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_HOME/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FFTW_HOME/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$I_MPI_ROOT/intel64/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$I_MPI_ROOT/intel64/lib/release/

cd /ptmp/mpa/hitesh/athena_plusplus_townsend/mixing_layer/

srun ./athena -i athinput_mixing.turb -t 05:55:00

echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo "Boom!"
