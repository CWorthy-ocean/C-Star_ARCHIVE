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

cstar_setupdir=$(pwd)

################################################################################
# MACHINE SPECIFIC COMMANDS AND VARIABLE SETTINGS
case "$1" in
    osx_arm64_gnu)
	# 1. Set up environment
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
	    
	    source $(conda info --base)/etc/profile.d/conda.sh
	    conda activate "${cstarenv}"
	    

	    # 2. CHECKOUT EXTERNALS
	    echo "getting external packages as specified in Externals.cfg"
	    ./manage_externals/checkout_externals

	    # Trim prefix in shell PS1 to just env dirname, not full path	    
	    conda config --set env_prompt '({name})'
	    # 3. SET ENVIRONMENT VARIABLES
	    conda env config vars set ROMS_ROOT="$(pwd)/externals/ucla-roms/" -p "${cstarenv}" > /dev/null
	    conda env config vars set MARBL_ROOT="$(pwd)/externals/MARBL/" -p "${cstarenv}" > /dev/null
	    conda env config vars set NETCDFHOME="${cstarenv}" -p "${cstarenv}" > /dev/null
	    conda env config vars set MPIHOME="${cstarenv}" -p "${cstarenv}" > /dev/null
	    conda activate "${cstarenv}" # The following are dependent on previous env vars so need to activate to see them
	    conda env config vars set LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NETCDFHOME/lib" -p "${cstarenv}" > /dev/null
	    conda env config vars set PATH="$ROMS_ROOT/Tools-Roms:$PATH" -p "${cstarenv}" > /dev/null
	    conda activate "${cstarenv}"
	    
	    
	fi
	# Alias command to go in .bashrc for future environment activation
 	CSTAR_ENV_ALIAS="conda activate ${cstarenv}"

	# 2. Set target for MARBL compilation
	MARBL_TARGET=gnu
	
	# 3. Distribute working ROMS makefiles for osx_arm64
	ls ${ROMS_ROOT}
	cd machine_specific_files/osx_arm64
	rsync -av ROMS_Makefiles/* ${ROMS_ROOT}
	;;
    
    sdsc_expanse_intel)
	if [ "$LMOD_SYSHOST" != "expanse" ];then
	    echo "You do not appear to be on the SDSC Expanse system.
	    If you believe you are reading this message in error, please raise an issue:
	    https://github.com/CWorthy-ocean/C-Star/issues/new
	    Exiting setup
	    "
	    exit 1
	else
	    echo "compiling on SDSC Expanse"
	fi
	#1. Set up environment

	# Set up within this script
	source machine_specific_files/sdsc_expanse/expanse_environment.sh
	#source ~/.bashrc
	ROMS_ROOT=$(pwd)/externals/ucla-roms/
	MARBL_ROOT=$(pwd)/externals/MARBL/

	# Alias command to go in bashrc for future environment activation
	CSTAR_ENV_ALIAS=\'"source $(pwd)/machine_specific_files/sdsc_expanse/expanse_environment.sh"\'

	# 2. Set target for MARBL compilation
	MARBL_TARGET=intel
	;;
    
    sdsc_expanse_gnu)
	echo "https://www.youtube.com/watch?v=ogwouE_Msd4"
	;;
   
esac

# UNIVERSAL COMMANDS USING ABOVE CONFIGURED ENVIRONMENTS AND VARIABLES
################################################################################

# 1. CHECKOUT EXTERNALS
echo "getting external packages as specified in Externals.cfg"
./manage_externals/checkout_externals

# 2. COMPILE MARBL
cd ${MARBL_ROOT}/src
sed -i -e "s#USEMPI=FALSE#USEMPI=TRUE#g" Makefile
echo "THE MARBL TARGET IS ${MARBL_TARGET}"

make "${MARBL_TARGET}"
cd ${cstar_setupdir}

# 3. COMPILE ROMS/NHMG and ROMS/TOOLS-ROMS LIBRARIES	
## i. make NHMG library
cd ${ROMS_ROOT}/Work
make nhmg 

## ii. make Tools-Roms
cd ${ROMS_ROOT}/Tools-Roms
make 
echo $ROMS_ROOT

# 4. Establish access to the C-Star environment in future
echo "All compilation steps successful!!" 
if ! grep -q "${CSTAR_ENV_ALIAS}" ~/.bashrc; then
    read -p \
	 "This setup script will now modify your .bashrc file to enable an environment for running C-Star.
      	 By default, an alias (cstar_env) will be added so that you can quickly set the environment up. Continue? (y/n)" \
	 continuestring
    
    case ${continuestring} in
	y|Y)
	    echo "#SET UP ENVIRONMENT FOR RUNNING C-STAR:" >> ~/.bashrc
	    echo "################################################################################" >> ~/.bashrc
	    echo "alias cstar_env=${CSTAR_ENV_ALIAS}" >> ~/.bashrc
	    echo "################################################################################" >> ~/.bashrc
	    
	    echo "Your ~/.bashrc file has been modified. In future, enter "
	    echo "cstar_env"
	    echo " to activate the environment for running C-Star."
	    echo "NOTE: This will not work until you either log back in or run"
	    echo "source ~/.bashrc"
	    ;;
	*)
	    echo "Your .bashrc file has not been modified. You will need to manually activate the C-Star environment using"
	    echo "${CSTAR_ENV_ALIAS}"
	    ;;
    esac
fi

