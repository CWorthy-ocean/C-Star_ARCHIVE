#!/bin/bash
set -e
#--------------------------------------------------------------------------------
# SCRIPT FOR FIRST TIME SETUP
# This program will acquire the necessary codebases, create an environment for
# your machine, and set up libraries necessary to run

# To set up a specific configuration, use the get_config script
#--------------------------------------------------------------------------------
if [ "$#" -eq 0 ]; then
    echo "please supply an argument and ask the dev to write a better error msg"
    exit 1
fi

# 1. CHECKOUT EXTERNALS
echo "getting external packages as specified in Externals.cfg"
./manage_externals/checkout_externals

cstar_setupdir=$(pwd)

case "$1" in
    osx_arm64_gnu)
	# 2. SETUP ENVIRONMENT
	# Check if env exists, install if not
	if [ -d conda_envs/cstar_gnu ] ;
	then echo "conda environment exists.
       	   	If you have already run this setup script and
       	       	 are experiencing problems with the environment
		 then remove it using"
	     echo "conda env remove -p $(pwd)/conda_envs/cstar_gnu"
	     echo "and try running this script again"
	else
	    cstarenv="$(pwd)/conda_envs/cstar_gnu"
	    conda env create -f conda_envs/cstar_gnu.yml --prefix="${cstarenv}"
	    echo "created conda environment"
	    source $(conda info --base)/etc/profile.d/conda.sh
	    conda activate "${cstarenv}"
	    echo "activated conda environment"
	    # Trim prefix in shell PS1 to just env dirname, not full path	    
	    conda config --set env_prompt '({name})'
	    # SET ENVIRONMENT VARIABLES
	    conda env config vars set ROMS_ROOT="$(pwd)/externals/ucla-roms/" -p "${cstarenv}" > /dev/null
	    conda env config vars set MARBL_ROOT="$(pwd)/externals/MARBL/" -p "${cstarenv}" > /dev/null
	    conda env config vars set NETCDFHOME="${cstarenv}" -p "${cstarenv}" > /dev/null
	    conda env config vars set MPIHOME="${cstarenv}" -p "${cstarenv}" > /dev/null
	    conda activate "${cstarenv}" # The following are dependent on previous env vars so need to activate to see them
	    conda env config vars set LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NETCDFHOME/lib" -p "${cstarenv}" > /dev/null
	    conda env config vars set PATH="$ROMS_ROOT/Tools-Roms:$PATH" -p "${cstarenv}" > /dev/null
	    conda activate "${cstarenv}"

	fi

	# 3. COMPILE MARBL
	cd ${MARBL_ROOT}/src
	gsed -i -e "s#USEMPI=FALSE#USEMPI=TRUE#g" Makefile
	make gnu
	cd ${cstar_setupdir}
	
	# 4. COMPILE ROMS/NHMG and ROMS/TOOLS-ROMS LIBRARIES
	# i. Distribute correct makefiles
	ls ${ROMS_ROOT}
	cd machine_specific_files/osx_arm64
	rsync -av ROMS_Makefiles/* ${ROMS_ROOT}
	
	# ii. make NHMG library
	cd ${ROMS_ROOT}/Work
	make nhmg 

	# iii. make Tools-Roms
	cd ${ROMS_ROOT}/Tools-Roms
	make # < FAILING HERE, DOESN'T SEEM TO BE LINKING NETCDF LIBRARY CORRECTLY
	
    ;;
    expanse_gnu)


	
    ;;
    expanse_intel)
    ;;



esac
