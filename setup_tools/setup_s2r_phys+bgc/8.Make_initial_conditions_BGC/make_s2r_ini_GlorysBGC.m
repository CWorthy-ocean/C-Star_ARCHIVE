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
    glorys_dir = '/glade/scratch/bachman/GLORYS/NA/2012/'
    glorys_mon_tr = [glorys_dir,'mercatorglorys12v1_gl12_mean_20120101.nc'];
    glorys_mon_vel= glorys_mon_tr;
    glorys_mon_ssh= glorys_mon_tr;
end
%%%%%
    % ROMS info
     romsdir    = '/glade/derecho/scratch/bachman/C-Star/setup_tools/setup_s2r_phys+bgc/1.Make_grid/';
     grdname    = [romsdir 'Iceland0_grd.nc'];
     ininame    = [romsdir 'Iceland0_ini.nc'];

%     %%%%% ref %%%%%
    pars.theta_s = 5.0;
    pars.theta_b = 2.0;
    pars.hc     = 300.0;
    pars.N      = 100;
    pars.scoord = 'new2012';    % child 'new' or 'old' type scoord

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    makebgc = 1 ; % turn on if bgc on
    %
    if makebgc==1
    %
    BGC_INI.data_woa18     = 1 ; BGC_INI.tracer_woa18     = [1 2 3 5] ;
    BGC_INI.data_CCSM      = 1 ; BGC_INI.tracer_CCSM      = [4] ;
    BGC_INI.data_SYnn      = 1 ; BGC_INI.tracer_SYnn      = [26 28] ;
    BGC_INI.data_glodapv2  = 1 ; BGC_INI.tracer_glodapv2  = [18 19] ;
    BGC_INI.data_SeaWiFS   = 1 ; BGC_INI.tracer_SeaWiFS   = [6] ;
    BGC_INI.data_Takahashi = 1 ; BGC_INI.tracer_Takahashi = [30] ;
    % if a data file is added , edit s2rhv_inbgc and create a new get_bgc_data.m file
    %                      /glade/scratch/bachman/ROMS_tools/code/BGC/forcings
    BGC_INI.file_bgc{1} = '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_n_annual_landfilled.nc' ;   % '/data/project1/data/WOA18/Nitrate/woa18_all_n*nc';
    BGC_INI.file_bgc{2} = '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_o_seasonal_landfilled.nc' ;  %'/data/project1/data/WOA18/Oxygen/woa18_all_o*nc';
    BGC_INI.file_bgc{3} = '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_p_annual_landfilled.nc' ; %'/data/project1/data/WOA18/Phosphate/woa18_all_p*nc';
    BGC_INI.file_bgc{4} = '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_i_seasonal_landfilled.nc'; %'/data/project1/data/WOA18/Silicate/woa18_all_i*nc';
    BGC_INI.file_bgc{5} = '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/ETHZ_L1_46.1k_Anth_Mon.pop.h.1958-2006_Fe_1981-2000_Feymonavg_rename_add_time_dimvar.nc';
    BGC_INI.file_bgc{6} = '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/n2ofromnn.nc' ;
    BGC_INI.file_bgc{7} = '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/GLODAPv2.2016b.TCO2.nc' ;
    BGC_INI.file_bgc{8} = '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/GLODAPv2.2016b.TAlk.nc' ;
    BGC_INI.file_bgc{9} = '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/SeaWiFS_CHL_MO_climatology_9km_landfill.nc';
    BGC_INI.file_bgc{10}= '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/n2o_satn2o_sat_mon_mmolpm3_fromWOA13_landfill.nc' ;
    BGC_INI.file_bgc{11}= '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/dpco2_climatology_T09_rev_landfill.nc';
    BGC_INI.file_bgc{12}= '/glade/derecho/scratch/bachman/ROMS_tools/DATASETS/BGC/pacmed_0p25_ini2000_5daysAVGsoda.nc';
    %
    BGC_INI.bgc_tracer = {'NO3','PO4','SiO3','Fe','O2','CHLA', ...
                  'SPC','SPCHL','SPFE','SPCACO3', ...
                  'DIATC','DIATCHL','DIATFE','DIATSI',...
                  'DIAZC','DIAZCHL','DIAZFE', ...
                  'DIC','Alk','ZOOC', ...
                  'DON','DONR','DOP','DOPR','DOFE'...
                  'N2O', 'N2O_SIDEN', 'N2O_ATM', 'N2O_NEV','pCO2','basindx'};

		 %{'NO3','PO4','SiO3','Fe','O2','CHLA', ...
                 % 'SPC','SPCHL','SPFE','SPCACO3', ...
                 % 'DIATC','DIATCHL','DIATFE','DIATSI',...
                 % 'DIAZC','DIAZCHL','DIAZFE', ...
                 % 'DIC_glodap','Alk_glodap','ZOOC', ...
                 % 'DON','DONR','DOP','DOPR','DOFE'...
                 % 'N2O', 'N2O_SIDEN', 'N2O_ATM', 'N2O_NEV','pCO2','basindx'}; %, ...
                 % 'DIC','Alk'}; % +PIC+CHL
    %  files to put in : 'i': ini/frc-file,
    %                    'v', value ,
    %                    'b' , build from other variables (_{Variable number})
    %                    's' , scaled from another variable (_{Variable number})

    bgc_tracer_list ;
    for trc=1:length(bgctracers_list.name)
      ind_trc = find(strcmp(bgctracers_list.name{trc}, BGC_INI.bgc_tracer));
      if isempty(ind_trc)
        BGC_INI.bgc_tracer{end+1} = bgctracers_list.name{trc};
      end
    end

    BGC_INI.bgc_frctype = {'i','i','i','i','i','i' ...
                  's_6','s_6','s_6','s_6', ...
                  's_6','s_6','s_6','s_6',...
                  's_6','s_6','s_6', ...
                  'i','i','s_6', ...
                  'v','v','v','v','v'...
                  'i', 'b_26', 'i', 'v', 'i', 'i', ...
                  'b_30','b_30'};

    for trc=length(BGC_INI.bgc_frctype)+1:length(bgctracers_list.name)
        BGC_INI.bgc_frctype{trc} = 'i';
    end

    % file to read for forcing if type 'i'; value to init  if 'v' or scale factor if 'b'
    BGC_INI.bgc_frcini = {BGC_INI.file_bgc{1},BGC_INI.file_bgc{3},BGC_INI.file_bgc{4},BGC_INI.file_bgc{5},BGC_INI.file_bgc{2},BGC_INI.file_bgc{9} ...
                  '3.375','0.675','1.35e-05','0.0675', ...
                  '0.2025','0.0675','1.35e-06','0.0675',...
                  '0.0375','0.0075','7.5e-07', ...
                  BGC_INI.file_bgc{7}, BGC_INI.file_bgc{8},'1.35', ...
                  '1','0.8','0.1','0.003','0.0001',...
                  BGC_INI.file_bgc{6}, ' ', BGC_INI.file_bgc{10}, '0.0001',BGC_INI.file_bgc{11},BGC_INI.file_bgc{12},...
                  ' ',' '};

    for trc=length(BGC_INI.bgc_frcini)+1:length(bgctracers_list.name)
        BGC_INI.bgc_frcini{trc} = ' ';
    end
    %
    %%%bgc_tracer_list ;
    %
    else
    BGC_INI = 0 ;
    bgctracers_list = 0 ;
    end
%
%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS END HERE
%---------------------------------------------------------------------------------------
%

%    if ~exist(ininame)
      disp(['Creating initial file: ' ininame]);
      r2r_create_ini(ininame,grdname,pars.N,pars,makebgc,BGC_INI,bgctracers_list);
%    end

     %if strcmp(SOURCE,'Soda')
     %s2r_hv_ini(soda_mon_tr,soda_mon_vel,soda_mon_ssh,grdname,ininame,pars)
     %elseif strcmp(SOURCE,'Glorys')
     %s2r_hv_ini(glorys_mon_tr,glorys_mon_vel,glorys_mon_ssh,grdname,ininame,pars)
     %end

    if makebgc==1
    s2r_hv_inibgc(BGC_INI,grdname,ininame,pars)
    end




