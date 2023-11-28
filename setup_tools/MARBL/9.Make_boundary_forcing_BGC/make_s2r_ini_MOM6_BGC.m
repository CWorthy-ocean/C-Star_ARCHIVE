%---------------------------------------------------------------------------------------
%
%  make_s2r_ini_MOM6_BGC
%
%  Generate boundary perimeter file from MOM6 data.
%
%  Scott Bachman, 2023 at [C]Worthy
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
%   File names:
%
%%%%% SOURCE = Soda or Glorys %%%%% Update s2r_hv_ini according`

   f_start = 1;
   f_end = 12;

   %  Data climatology file names:
   src_dir = '/glade/derecho/scratch/mlevy/archive/g.e23b15.TL319_t061.G1850MOMMARBL_JRA.012/ocn/hist/';
   src = [src_dir,'g.e23b15.TL319_t061.G1850MOMMARBL_JRA.012.mom6.h.bgc.*.nc'];
   src_files  = dir(src) ;
   src_files = src_files(f_start:f_end);


    % ROMS info
     romsdir    = '/glade/cheyenne/scratch/dafydd/ucla-roms/run/roms_marbl_BATS/ROMS_tools/1.Make_grid/';
     outdir = '/glade/cheyenne/scratch/bachman/C-Star/setup_tools/MARBL/9.Make_boundary_forcing_BGC/';
     grdname    = [romsdir 'roms_grd.nc'];
     bryname    = [outdir 'roms_bryBGC.nc'];

    %%%%% ref %%%%%
    pars.theta_s = 5.0;
    pars.theta_b = 2.0;
    pars.hc     = 300.0;
    pars.N      = 100;
    pars.scoord = 'new2012';    % child 'new' or 'old' type scoord

    obcflag        = [1 1 1 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    BGC_INI.bgc_tracer{1} = 'PO4';
    BGC_INI.bgc_tracer{end+1} = 'NO3';
    BGC_INI.bgc_tracer{end+1} = 'SiO3';
    BGC_INI.bgc_tracer{end+1} = 'NH4';
    BGC_INI.bgc_tracer{end+1} = 'Fe';
    BGC_INI.bgc_tracer{end+1} = 'Lig';
    BGC_INI.bgc_tracer{end+1} = 'O2';
    BGC_INI.bgc_tracer{end+1} = 'DIC';
    BGC_INI.bgc_tracer{end+1} = 'DIC_ALT_CO2';
    BGC_INI.bgc_tracer{end+1} = 'ALK';
    BGC_INI.bgc_tracer{end+1} = 'ALK_ALT_CO2';
    BGC_INI.bgc_tracer{end+1} = 'DOC';
    BGC_INI.bgc_tracer{end+1} = 'DON';
    BGC_INI.bgc_tracer{end+1} = 'DOP';
    BGC_INI.bgc_tracer{end+1} = 'DOPr';
    BGC_INI.bgc_tracer{end+1} = 'DONr';
    BGC_INI.bgc_tracer{end+1} = 'DOCr';
    BGC_INI.bgc_tracer{end+1} = 'zooC';
    BGC_INI.bgc_tracer{end+1} = 'spChl';
    BGC_INI.bgc_tracer{end+1} = 'spC';
    BGC_INI.bgc_tracer{end+1} = 'spP';
    BGC_INI.bgc_tracer{end+1} = 'spFe';
    BGC_INI.bgc_tracer{end+1} = 'spCaCO3';
    BGC_INI.bgc_tracer{end+1} = 'diatChl';
    BGC_INI.bgc_tracer{end+1} = 'diatC';
    BGC_INI.bgc_tracer{end+1} = 'diatP';
    BGC_INI.bgc_tracer{end+1} = 'diatFe';
    BGC_INI.bgc_tracer{end+1} = 'diatSi';
    BGC_INI.bgc_tracer{end+1} = 'diazChl';
    BGC_INI.bgc_tracer{end+1} = 'diazC';
    BGC_INI.bgc_tracer{end+1} = 'diazP';
    BGC_INI.bgc_tracer{end+1} = 'diazFe';

    bgc_tracer_list ;
    for trc=1:length(bgctracers_list.name)
      ind_trc = find(strcmp(bgctracers_list.name{trc}, BGC_INI.bgc_tracer));
      if isempty(ind_trc)
        disp('WARNING: Tracer ', bgctracers_list.name{trc}, ' in tracers list but will not be initialized.')
      end
    end

    for t = 1:(f_end-f_start+1)
      BRYtime.source{t} = [src_dir src_files(t).name];
    end

%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS END HERE
%---------------------------------------------------------------------------------------

    disp(['Creating boundary file: bryname' bryname]);
    create_bry(bryname,grdname,obcflag,pars,BRYtime,BGC_INI,bgctracers_list);

    s2r_hvbgc(BGC_INI,grdname,bryname,pars,obcflag,BRYtime)




