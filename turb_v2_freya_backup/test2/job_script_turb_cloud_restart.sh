#!/bin/sh
#
#SBATCH -J turb_rs
#SBATCH -o turb_rs."%j".out
#SBATCH -e turb_rs."%j".err
#SBATCH --mail-user hitesh@mpa-garching.mpg.de
#SBATCH --mail-type=ALL
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time="00:29:00"

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

cd /ptmp/hitesh/athena_fork_turb_box/turb_v2/test2/

srun ./athena_turb -r Turbulence/Turb.00010.rst -i athinput.turb_v2_cloud -t 00:28:00 


echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo "Boom!"
