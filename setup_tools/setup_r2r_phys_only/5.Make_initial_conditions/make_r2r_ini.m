%---------------------------------------------------------------------------------------
%
%  make_s2r
%
%  Generate boundary perimeter file from WOA and SSH  data.
%
%  Note that when run this script it tests for the presence of a .mat file
%  which contains various interpolation coefficients related to your child
%  and parent grids.  If the .mat file is not there it will calculate the coefficients
%
%
%  Jeroen Molemaker and Evan Mason in 2007-2009 at UCLA
%
%---------------------------------------------------------------------------------------
clear all
close all
disp(' ')
%---------------------------------------------------------------------------------------
%  USER-DEFINED VARIABLES & OPTIONS START HERE
%---------------------------------------------------------------------------------------
%
%---------------------------------------------------------------------------------------
%  1.  GENERAL
%---------------------------------------------------------------------------------------
%
%   Parent...
%
     parscd.file    = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/CT1_2019-2020/HOLD_RST/CT1_rst.20190301120000.nc';
     pargrd = '/glade/derecho/scratch/bachman/UCLA-ROMS/Work/CT/CT1/INPUT/2019-2020/CT1_grd.nc' ;
     parscd.N       = 100 ;
     parscd.theta_s = 5.0;
     parscd.theta_b = 2.0;
     parscd.hc      = 300 ;
     parscd.tind    = 2;            % frame number in parent file
     parscd.scoord = 'new2012';    % child 'new' or 'old' type scoord

%%%%% child
    romsdir    = '/glade/cheyenne/scratch/bachman/C-Star/setup_tools/setup_r2r_phys_only/1.Make_grid/';
    chdgrd    = [romsdir 'CT2_grd.nc'];
    chdini    = [romsdir 'CT2_ini.nc'];
    chdscd.theta_s = 5.0;
    chdscd.theta_b = 2.0;
    chdscd.hc     = 300.0;
    chdscd.N      = 100;
    chdscd.scoord = 'new2012';    % child 'new' or 'old' type scoord

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS END HERE
%---------------------------------------------------------------------------------------
%

    if ~exist(chdini)
      dateref = datenum(2000,1,1);
      %ini_time = ( datenum(1950, 1, 1) + ncread(parscd.file,'ocean_time')/24 - dateref ) %* 24 * 3600
      restart_time = ncread(parscd.file,'ocean_time');
      ini_time = restart_time(parscd.tind);
      %ini_time = ncread(parscd.file,'ocean_time');
      disp(['Creating initial file: ' chdini]);
      %r2r_create_ini(chdini,chdgrd,chdscd.N,chdscd)
      r2r_create_ini(chdini,chdgrd,chdscd.N,chdscd, ini_time)
    end

    r2r_make_ini(pargrd, parscd.file, chdgrd, chdini, chdscd,parscd,parscd.scoord,chdscd.scoord)



