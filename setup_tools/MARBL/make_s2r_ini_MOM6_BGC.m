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
    % ROMS info
    romsdir    = '/glade/derecho/scratch/bachman/C-Star/setup_tools/setup_s2r_phys+bgc/1.Make_grid/';
    grdname    = [romsdir 'Iceland0_grd.nc'];
    ininame    = [romsdir 'Iceland0_ini_bgc.nc'];

    %%%%% ref %%%%%
    pars.theta_s = 5.0;
    pars.theta_b = 2.0;
    pars.hc     = 300.0;
    pars.N      = 100;
    pars.scoord = 'new2012';    % child 'new' or 'old' type scoord

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    BGC_INI.bgc_tracer{1} = 'PO4';
    BGC_INI.bgc_tracer{2} = 'NO3';
    BGC_INI.bgc_tracer{3} = 'SiO3';
    BGC_INI.bgc_tracer{4} = 'NH4';
    BGC_INI.bgc_tracer{5} = 'Fe';
    BGC_INI.bgc_tracer{6} = 'Lig';
    BGC_INI.bgc_tracer{7} = 'O2';
    BGC_INI.bgc_tracer{8} = 'DIC';
    BGC_INI.bgc_tracer{9} = 'DIC_ALT_CO2';
    BGC_INI.bgc_tracer{10} = 'ALK';
    BGC_INI.bgc_tracer{11} = 'ALK_ALT_CO2';
    BGC_INI.bgc_tracer{12} = 'DOC';
    BGC_INI.bgc_tracer{13} = 'DON';
    BGC_INI.bgc_tracer{14} = 'DOP';
    BGC_INI.bgc_tracer{15} = 'DOPr';
    BGC_INI.bgc_tracer{16} = 'DONr';
    BGC_INI.bgc_tracer{17} = 'DOCr';
    BGC_INI.bgc_tracer{18} = 'zooC';
    BGC_INI.bgc_tracer{19} = 'spChl';
    BGC_INI.bgc_tracer{20} = 'spC';
    BGC_INI.bgc_tracer{21} = 'spP';
    BGC_INI.bgc_tracer{22} = 'spFe';
    BGC_INI.bgc_tracer{23} = 'spCaCO3';
    BGC_INI.bgc_tracer{24} = 'diatChl';
    BGC_INI.bgc_tracer{25} = 'diatC';
    BGC_INI.bgc_tracer{26} = 'diatP';
    BGC_INI.bgc_tracer{27} = 'diatFe';
    BGC_INI.bgc_tracer{28} = 'diatSi';
    BGC_INI.bgc_tracer{29} = 'diazChl';
    BGC_INI.bgc_tracer{30} = 'diazC';
    BGC_INI.bgc_tracer{31} = 'diazP';
    BGC_INI.bgc_tracer{32} = 'diazFe';

    BGC_INI.source{1} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{2} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{3} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{4} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{5} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{6} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{7} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{8} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{9} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{10} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{11} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{12} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{13} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{14} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{15} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{16} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{17} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{18} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{19} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{20} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{21} = '/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';
    BGC_INI.source{22} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{23} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{24} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{25} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{26} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{27} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{28} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{29} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{30} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{31} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';
    BGC_INI.source{32} = '/glade/work/mlevy/cesm_inputdata/MOM_IC_1.nc';

    bgc_tracer_list ;
    for trc=1:length(bgctracers_list.name)
      ind_trc = find(strcmp(bgctracers_list.name{trc}, BGC_INI.bgc_tracer));
      if isempty(ind_trc)
        disp('WARNING: Tracer ', bgctracers_list.name{trc}, ' in tracers list but will not be initialized.')
      end
    end

%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS END HERE
%---------------------------------------------------------------------------------------

    disp(['Creating initial file: ' ininame]);
    create_ini(ininame,grdname,pars.N,pars,BGC_INI,bgctracers_list);

    s2r_hv_inibgc(BGC_INI,grdname,ininame,pars)




