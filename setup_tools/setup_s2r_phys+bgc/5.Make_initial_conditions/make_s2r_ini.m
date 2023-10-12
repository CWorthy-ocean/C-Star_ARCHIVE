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
%   Data climatology file names:
%

%%%%% SOURCE = Soda or Glorys %%%%% Update s2r_hv_ini according`
SOURCE = 'Glorys' ;

if strcmp(SOURCE,'Soda')
    hbl_dir = '/data/project2/data/SODA/soda_3.3.1/5day_orig/';
    hbl_mon_data = [hbl_dir,'soda3.3.1_5dy_ocean_or_2000_01_03.nc'];
    %
    ssh_dir = '/data/project2/data/SODA/soda_3.3.1/5day_orig/';
    ssh_mon_data  = [ssh_dir,'soda3.3.1_5dy_ocean_or_2000_01_03.nc'];
    %
    soda_dir = '/data/project2/data/SODA/soda_3.3.1/5day_orig/';
    soda_mon_tr = [soda_dir,'soda3.3.1_5dy_ocean_or_2000_01_03.nc'];
    soda_mon_vel= [soda_dir,'soda3.3.1_5dy_ocean_or_2000_01_03.nc'];
    soda_mon_ssh= [soda_dir,'soda3.3.1_5dy_ocean_or_2000_01_03.nc'];

elseif strcmp(SOURCE,'Glorys')

    glorys_dir = ['/glade/scratch/bachman/GLORYS/NA/2012/']
    soda_mon_tr  = [glorys_dir,'mercatorglorys12v1_gl12_mean_20120102.nc'];
    soda_mon_vel = soda_mon_tr;
    soda_mon_ssh = soda_mon_tr;

end

%%%%%
    % ROMS info

    grdname = ['/glade/scratch/bachman/ROMS_tools/Iceland0_BGC2/1.Make_grid/Iceland0_grd.nc'];
    ininame = ['/glade/scratch/bachman/ROMS_tools/Iceland0_BGC2/1.Make_grid/Iceland0_ini.nc'];

    pars.theta_s = 5.0;
    pars.theta_b = 2.0;
    pars.hc     = 300.0;
    pars.N      = 100;
    pars.scoord = 'new2012';    % child 'new' or 'old' type scoord

    dateref = datenum(2000,1,1)
    % for GLORYS:
    ini_time = ( datenum(1950, 1, 1) + ncread(soda_mon_tr,'time')/24 - dateref ) * 24 * 3600 ;

%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS END HERE
%---------------------------------------------------------------------------------------

    if ~exist(ininame)
      disp(['Creating initial file: ' ininame]);
      r2r_create_ini(ininame,grdname,pars.N,pars,ini_time);
    end

    s2r_hv_ini(soda_mon_tr,soda_mon_vel,soda_mon_ssh,grdname,ininame,pars,ini_time,SOURCE)


