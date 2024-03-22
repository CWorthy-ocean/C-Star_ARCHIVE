# .bashrc
if [ -n "$PS1" ];then 
    # Source global definitions
    if [ -f /etc/bashrc ]; then
	. /etc/bashrc
    fi
    
    # User specific environment
    if ! [[ "$PATH" =~ "$HOME/.local/bin:$HOME/bin:" ]]
    then
	PATH="$HOME/.local/bin:$HOME/bin:$PATH"
    fi
    export PATH
    export PATH="/cm/shared/apps/sdsc/galyleo:${PATH}"
    
    # Uncomment the following line if you don't like systemctl's auto-paging feature:
    # export SYSTEMD_PAGER=
    
    # User specific aliases and functions
    
    #module load gcc
    module load slurm
    
    # Account Managing
    module load sdsc
    alias budget='expanse-client project edf100 -p'
    #lfs quota -g edf100 -h /expanse/lustre/projects/edf100
    
    # additions:
    
    # - expanse user-defined path to where you have cloned the repo code:
    export ROMS_ROOT='/expanse/lustre/projects/edf100/dstephenson/my-ucla-roms/'
    
    export Q='/cm/shared/apps/sdsc/current/bin/expanse-client user'
    export SACCT_FORMAT="jobid,jobname,partition,account,alloccpus,state,elapsed,timelimit"
    
    # Expanse non-user defined (same for everyone)
    # - Prerequisite modules for netcdf (as told by using 'module spider netcdf')
    module load cpu/0.15.4  intel/19.1.1.217  mvapich2/2.3.4
    # - netCDF-c
    module load netcdf-c/4.7.4
    # - netCDF-Fortran
    module load netcdf-fortran/4.5.3
    # - ncview
    module load ncview/2.1.8
    
    # - adds roms tools to path (you still need to build them first)
    export PATH=$PATH:$ROMS_ROOT/Tools-Roms
    
    # - set roms' environment variables to match expanse module paths:
    # 1 - ucla-roms compilation
    export NETCDFHOME=$NETCDF_FORTRANHOME
    export MPIHOME=$MVAPICH2HOME
    # 2 - older roms verions compilation
    export NETCDF=$NETCDF_FORTRANHOME
    export MPI_ROOT=$MVAPICH2HOME
    # 3 - MARBL
    export MARBL_ROOT=/home/dstephenson/MARBL
    # - Show my current submitted jobs
    alias SQ='squeue -u dstephenson'
    alias ll='ls -l'
    alias rmtil="rm *\~"
    alias rmhas="rm \#*"
    alias tailout="tail -f $(ls -t *.out | head -1)"
    alias viewout="emacs  $(ls -t *.out | head -1)"
    alias cdl='cd /expanse/lustre/projects/edf100/dstephenson'
    export SCRATCH='/expanse/lustre/scratch/dstephenson/temp_project'
    export RUN='/expanse/lustre/projects/edf100/dstephenson'
    ### Obsolete Comet aliases ###
    #alias matlabnd='module load matlab; matlab -nodesktop -nodisplay -nosplash'
    #alias intseshD='srun --partition=debug --pty --nodes=1 --ntasks-per-node=2 -t 00:30:00 --wait=0 --export=ALL /bin/bash'
    #alias intseshC='srun --partition=compute --pty --nodes=1 --ntasks-per-node=1 -t 04:00:00 --wait=0 --export=ALL /bin/bash'
    #alias intsesh40='srun --partition=compute --pty --nodes=1 --mem=50G --ntasks-per-node=24 -t 24:00:00 --wait=0 --export=ALL /bin/bash'
    #alias intseshlarge='srun --partition=large-shared --mem=100G --pty --nodes=1 --ntasks-per-node=12 -t 24:00:00 --export=ALL /bin/bash'
    
    conda_init=bash
    
    # >>> conda initialize >>>
    # !! Contents within this block are managed by 'conda init' !!
    __conda_setup="$('/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen/gcc-8.3.1/anaconda3-2020.11-da3i7hmt6bdqbmuzq6pyt7kbm47wyrjp/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
	eval "$__conda_setup"
    else
	if [ -f "/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen/gcc-8.3.1/anaconda3-2020.11-da3i7hmt6bdqbmuzq6pyt7kbm47wyrjp/etc/profile.d/conda.sh" ]; then
	    # . "/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen/gcc-8.3.1/anaconda3-2020.11-da3i7hmt6bdqbmuzq6pyt7kbm47wyrjp/etc/profile.d/conda.sh"  # commented out by conda initialize
	    echo ""
	else
	    # export PATH="/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen/gcc-8.3.1/anaconda3-2020.11-da3i7hmt6bdqbmuzq6pyt7kbm47wyrjp/bin:$PATH"  # commented out by conda initialize
	    echo ""
	fi
    fi
    unset __conda_setup
    # <<< conda initialize <<<
    
    conda init bash
    
    
    if [[ ${SSH_TTY} == *"/dev/pts/"* ]];then
	green="\001\$(tput setaf 47)\002" #replaces \[ with \001 and \] with \002                                                                                                                                          
	blue="\001$(tput setaf 51)\002"
	red="\001$(tput setaf 196)\002"
	dim="\001$(tput dim)\002"
	reset="\001$(tput sgr0)\002"
	#    if $(echo ${HOSTNAME} | grep -q "cheyenne"); then
	PS1="$dim\t |" # time in hh:mm:ss
	PS1+="$blue \W ${reset}|" #working directory
	PS1+="$blue\u\$${reset}"
	#    elif $(echo ${HOSTNAME} | grep -q "casper"); then
	#        PS1="$dim\t |" # time in hh:mm:ss                                                                                                                                                                            
	#        PS1+="$red \W ${reset}|" #working directory                                                                                                                                                                  
	#        PS1+="$red\u\$${reset}"
	#    fi
	export PS1
	unset green blue dim reset
    fi
fi
