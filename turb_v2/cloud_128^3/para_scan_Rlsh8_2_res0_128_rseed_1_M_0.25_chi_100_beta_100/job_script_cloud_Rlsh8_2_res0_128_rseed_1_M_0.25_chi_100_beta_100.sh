#!/bin/sh
#
#SBATCH -J mhd80_1
#SBATCH -o mhd80_1."%j".out
#SBATCH -e mhd80_1."%j".err
#SBATCH --mail-user hitesh@mpa-garching.mpg.de
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:59:00

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

cd /ptmp/hitesh/athena_fork_turb_box/turb_v2/para_scan_Rlsh8_2_res0_128_rseed_1_M_0.25_chi_100_beta_100/

srun ./athena_turb -r ../Turbulence/para_scan_Rlsh8_2_res0_128_rseed_1_M_0.25_beta_100/Turb.final.rst -i athinput_cloud_Rlsh8_2_res0_128_rseed_1_M_0.25_chi_100_beta_100.turb -t 01:55:00

echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo "Boom!"
