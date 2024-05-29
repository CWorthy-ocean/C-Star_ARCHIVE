module restore

# Load C-star related modules
module load cpu
module load cray-hdf5/1.12.2.9
module load cray-netcdf/4.9.0.9

export MPIHOME=/opt/cray/pe/mpich/8.1.28/ofi/gnu/12.3/
export NETCDFHOME=/opt/cray/pe/netcdf/4.9.0.9/gnu/12.3/

export PATH=$PATH:$ROMS_ROOT/Tools-Roms
export PATH=$PATH:$NETCDFHOME/bin

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDFHOME/lib
export LIBRARY_PATH=$LIBRARY_PATH:$NETCDFHOME/lib
