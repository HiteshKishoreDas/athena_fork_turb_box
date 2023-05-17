#!/bin/sh
#
#SBATCH -J hyd00_1
#SBATCH -o hyd00_1."%j".out
#SBATCH -e hyd00_1."%j".err
#SBATCH --mail-user hitesh@mpa-garching.mpg.de
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=11:49:00

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

cd /ptmp/hitesh/athena_fork_turb_box/turb_v2/para_scan_Rlsh0_100_res0_128_rseed_1_M_0.5_hydro/

srun ./athena_turb -i athinput_turb_Rlsh0_100_res0_128_rseed_1_M_0.5_hydro.turb -t 11:35:00

echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo "Boom!"
