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

for yyyy=2012:2012

 clear list_soda_hbl
 clear list_soda_tr
 clear list_soda_vel
 clear list_soda_ssh
 clear BRYtime

%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS START HERE
%---------------------------------------------------------------------------------------
%
%---------------------------------------------------------------------------------------
% 1.  GENERAL
%---------------------------------------------------------------------------------------

   day_start = 1;
   day_end = 365;

   %  Data climatology file names:
   glorys_dir = '/glade/scratch/bachman/GLORYS/NA/2012/';
   glorys_mon_tr = [glorys_dir,'mercatorglorys12v1_gl12_mean_2012*.nc'];
   glorys_mon_vel= glorys_mon_tr;
   glorys_mon_ssh= glorys_mon_tr;

    % ROMS info
     romsdir    = '/glade/scratch/bachman/ROMS_tools/Iceland0_BGC2/1.Make_grid/';
     grdname    = [romsdir 'Iceland0_grd.nc'];
     bryname    = [romsdir 'Iceland0_bry_2012.nc'];

%     %%%%% ref %%%%%
    pars.theta_s = 5.0;
    pars.theta_b = 2.0;
    pars.hc     = 300.0;
    pars.N      = 100;
    pars.scoord = 'new2012';    % child 'new' or 'old' type scoord

    % OBC flags and time
    obcflag        = [1 1 1 1];    % open boundaries flag (1=open , [S E N W])
    bry_cycle      = 365.25;       % 0 means no cycle
    mdays = [31 28 31 30 31 30 31 31 30 31 30 31] ;
    DT             = 1 ;           % time step of the forcing file (days)
    extraband      = 1  ;          % extra value at the boundary
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
    %
    BGC_INI.file_bgc{1} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_n_annual_landfilled.nc' ;   % '/data/project1/data/WOA18/Nitrate/woa18_all_n*nc';
    BGC_INI.file_bgc{2} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_o_seasonal_landfilled.nc' ;  %'/data/project1/data/WOA18/Oxygen/woa18_all_o*nc';
    BGC_INI.file_bgc{3} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_p_annual_landfilled.nc' ; %'/data/project1/data/WOA18/Phosphate/woa18_all_p*nc';
    BGC_INI.file_bgc{4} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_i_seasonal_landfilled.nc'; %'/data/project1/data/WOA18/Silicate/woa18_all_i*nc';
    BGC_INI.file_bgc{5} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/ETHZ_L1_46.1k_Anth_Mon.pop.h.1958-2006_Fe_1981-2000_Feymonavg_rename_add_time_dimvar.nc';
    BGC_INI.file_bgc{6} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/n2ofromnn.nc' ;
    BGC_INI.file_bgc{7} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/GLODAPv2.2016b.TCO2.nc' ;
    BGC_INI.file_bgc{8} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/GLODAPv2.2016b.TAlk.nc' ;
    BGC_INI.file_bgc{9} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/SeaWiFS_CHL_MO_climatology_9km_landfill.nc';
    BGC_INI.file_bgc{10}= '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/n2o_satn2o_sat_mon_mmolpm3_fromWOA13_landfill.nc' ;
    BGC_INI.file_bgc{11}= '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/dpco2_climatology_T09_rev_landfill.nc';
    BGC_INI.file_bgc{12}= '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/pacmed_0p25_zclm2.nc';
    %
    BGC_INI.bgc_tracer = {'NO3','PO4','SiO3','Fe','O2','CHLA', ...
                  'SPC','SPCHL','SPFE','SPCACO3', ...
                  'DIATC','DIATCHL','DIATFE','DIATSI',...
                  'DIAZC','DIAZCHL','DIAZFE', ...
                  'DIC_glodap','Alk_glodap','ZOOC', ...
                  'DON','DONR','DOP','DOPR','DOFE'...
                  'N2O', 'N2O_SIDEN', 'N2O_ATM', 'N2O_NEV','pCO2','basindx', ...
                  'DIC','Alk'}; % +PIC+CHL

    bgc_tracer_list;
    for trc=1:length(bgctracers_list.name)
      ind_trc = find(strcmp(bgctracers_list.name{trc}, BGC_INI.bgc_tracer));
      if isempty(ind_trc)
        BGC_INI.bgc_tracer{end+1} = bgctracers_list.name{trc};
      end
    end

    %  files to put in : 'i': ini/frc-file,
    %                    'v', value ,
    %                    'b' , build from other variables (_{Variable number})
    %                    's' , scaled from another variable (_{Variable number})
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
    rep_tr  = dir(glorys_mon_tr ) ;
    rep_vel = dir(glorys_mon_vel) ;
    rep_ssh = dir(glorys_mon_ssh) ;

    rep_tr = rep_tr(day_start:day_end);
    rep_vel = rep_vel(day_start:day_end);
    rep_ssh = rep_ssh(day_start:day_end);

    for l=1:length(rep_tr)
    list_soda_tr {l} = [glorys_dir rep_tr(l,1).name ];
    list_soda_vel{l} = [glorys_dir rep_vel(l,1).name];
    list_soda_ssh{l} = [glorys_dir rep_ssh(l,1).name];
    end
    if (extraband>0)
        if extraband>1
            disp('ERROR extraband>1 not coded') ; error ;
        end
        if (yyyy==1994)
            disp('ERROR extraband init to reWork for 1994') ; error ;
            list_soda_ssh = ['/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_1980_01_03.nc' ...
                            list_soda_ssh ...
                            '/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_1981_01_02.nc'];
            list_soda_tr  = ['/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_1980_01_03.nc' ...
                            list_soda_tr ...
                            '/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_1981_01_02.nc'];
            list_soda_vel = ['/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_1980_01_03.nc' ...
                            list_soda_vel ...
                            '/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_1981_01_02.nc'];
        elseif (yyyy==2016)
            disp('ERROR extraband init to reWork for 2016') ; error ;
            list_soda_ssh = ['/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_2014_12_30.nc' ...
                            list_soda_ssh ...
                            '/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_2015_12_30.nc'];
            list_soda_tr  = ['/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_2014_12_30.nc' ...
                            list_soda_tr ...
                            '/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_2015_12_30.nc'];
            list_soda_vel = ['/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_2014_12_30.nc' ...
                            list_soda_vel ...
                            '/data/project2/data/SODA/soda_3.3.1/5day_orig/soda3.3.1_5dy_ocean_or_2015_12_30.nc'];
        else
            %%% Need to uncomment/edit this if we want more than one year at a time
            %{
            file_prev = dir([glorys_dir,'glorys12v1_Y' num2str(yyyy-1) '*.nc']) ;
            file_prev = [glorys_dir  file_prev(end,1).name];
            file_next = dir([glorys_dir,'glorys12v1_Y' num2str(yyyy+1) '*.nc']) ;
            file_next = [glorys_dir  file_next(1,1).name];
            list_soda_ssh = [file_prev list_soda_ssh file_next];

            file_prev = dir([glorys_dir,'glorys12v1_Y' num2str(yyyy-1) '*.nc']) ;
            file_prev = [glorys_dir  file_prev(end,1).name];
            file_next = dir([glorys_dir,'glorys12v1_Y' num2str(yyyy+1) '*.nc']) ;
            file_next = [glorys_dir  file_next(1,1).name];
            list_soda_tr = [file_prev list_soda_tr file_next];

            file_prev = dir([glorys_dir,'glorys12v1_Y' num2str(yyyy-1) '*.nc']) ;
            file_prev = [glorys_dir  file_prev(end,1).name];
            file_next = dir([glorys_dir,'glorys12v1_Y' num2str(yyyy+1) '*.nc']) ;
            file_next = [glorys_dir  file_next(1,1).name];
            list_soda_vel = [file_prev list_soda_vel file_next];
            %}
        end
    end

    %%% The character indices of month and day need to be hardcoded here
    %%% BRY time
    BRYtime.mdays = [31 28 31 30 31 30 31 31 30 31 30 31] ;
    file1 = list_soda_tr{extraband+1};
    t1 = str2num(file1(74:75)) ;
    for t=1:length(list_soda_tr)
        BRYtime.time(t) = t1 + DT*(t-1.5-extraband) ;
        file1 = list_soda_tr{t};
        BRYtime.month(t) = str2num(file1(72:73)) ;
        BRYtime.day(t)   = str2num(file1(74:75)) ;
        if BRYtime.day(t)<=BRYtime.mdays(BRYtime.month(t))/2
        BRYtime.interp_ratio(t) = BRYtime.day(t)/BRYtime.mdays(BRYtime.month(t)) + 0.5 ;
        else
        BRYtime.interp_ratio(t) = 1.5 - BRYtime.day(t)/BRYtime.mdays(BRYtime.month(t)) ;
        end
    end
    BRYtime.tstart = BRYtime.time(extraband+1);
    BRYtime.tend   = BRYtime.time(length(BRYtime.time)-extraband);
    BRYtime.cycle = bry_cycle ;

    % Create the bry file
    %if ~exist(bryname)
      disp(['Creating boundary file: bryname' bryname]);
      r2r_create_bry(bryname,grdname,obcflag,pars,BRYtime,makebgc,BGC_INI,bgctracers_list);
    %end

    for days =1:(day_end - day_start + 1)  %%%length(list_soda_tr)

    disp(['%%----> Working on days ' num2str(days) '/' num2str(length(list_soda_tr)) ' yy = ' num2str(yyyy) ' <----%%'])
    %s2r_hv_glorys(list_soda_tr{days},list_soda_vel{days},list_soda_ssh{days},grdname,bryname,days,pars,obcflag);
    if makebgc==1
       s2r_hvbgc(BGC_INI,grdname,bryname,days,pars,obcflag,BRYtime);
    end

    end

    if makebgc==1
       s2r_hvbgc_recalc(BGC_INI,grdname,bryname,pars,obcflag);
    end


end
