module purge

module load PrgEnv-intel
module load cray-hdf5
module load cray-netcdf

export MPIHOME=/opt/cray/pe/mpich/8.1.28/ofi/intel/2022.1/
export NETCDFHOME=/opt/cray/pe/netcdf/4.9.0.9/intel/2023.2/

export PATH=$PATH:$NETCDFHOME/bin

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDFHOME/lib
export LIBRARY_PATH=$LIBRARY_PATH:$NETCDFHOME/lib
