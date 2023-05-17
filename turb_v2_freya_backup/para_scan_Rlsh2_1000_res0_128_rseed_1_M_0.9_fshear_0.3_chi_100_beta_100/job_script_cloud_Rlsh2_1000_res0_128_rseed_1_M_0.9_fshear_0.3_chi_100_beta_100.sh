#!/bin/sh
#
#SBATCH -J m20M21
#SBATCH -o m20M21."%j".out
#SBATCH -e m20M21."%j".err
#SBATCH --mail-user hitesh@mpa-garching.mpg.de
#SBATCH --partition=p.24h
#SBATCH --mail-type=ALL
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --time=23:59:00

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

cd /ptmp/mpa/hitesh/athena_fork_turb_box/turb_v2/para_scan_Rlsh2_1000_res0_128_rseed_1_M_0.9_fshear_0.3_chi_100_beta_100/

cp ../Turbulence/para_scan_Rlsh2_1000_res0_128_rseed_1_M_0.9_fshear_0.3_beta_100/Turb.hst ./
# mv ../para_scan_Rlsh2_1000_res0_128_rseed_1_M_0.9_fshear_0.3_beta_100 ../Turbulence/

srun ./athena_turb -r ../Turbulence/para_scan_Rlsh2_1000_res0_128_rseed_1_M_0.9_fshear_0.3_beta_100/Turb.final.rst -i athinput_cloud_Rlsh2_1000_res0_128_rseed_1_M_0.9_fshear_0.3_chi_100_beta_100.turb -t 23:55:00

echo "Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo "Boom!"
