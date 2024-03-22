#!/bin/bash

################################################################################
# MODULES
module purge
# - Prerequisite modules for netcdf (as told by using 'module spider netcdf')
module load cpu/0.15.4  intel/19.1.1.217  mvapich2/2.3.4
# - netCDF-c
module load netcdf-c/4.7.4
# - netCDF-Fortran
module load netcdf-fortran/4.5.3
# - ncview
module load ncview/2.1.8
################################################################################
# ENVIRONMENT
export ROMS_ROOT=$(dirname $(pwd))
export MARBL_ROOT=$ROMS_ROOT/manage_externals/MARBL

# - set roms' environment variables to match expanse module paths:
# 1 - ucla-roms compilation
export NETCDFHOME=$NETCDF_FORTRANHOME
export MPIHOME=$MVAPICH2HOME
# 2 - older roms verions compilation
export NETCDF=$NETCDF_FORTRANHOME
export MPI_ROOT=$MVAPICH2HOME

################################################################################
