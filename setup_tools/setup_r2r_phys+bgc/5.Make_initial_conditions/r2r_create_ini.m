function r2r_create_ini(inifile,gridfile,N,chdscd,ini_time)
%
%   Input:
%
%   inifile      Netcdf initial file name (character string).
%   gridfile     Netcdf grid file name (character string).
%   clobber      Switch to allow writing over an existing
%                file (character string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Read the grid file
%
h       = ncread(gridfile,'h')';
[Mp,Lp] = size(h);
L       = Lp - 1 ;
M       = Mp - 1 ;
Np      = N  + 1 ;

%
%
%  Create variables and attributes
%
nccreate(inifile,'tstart','Dimensions',{'one',1},'datatype','single');
ncwriteatt(inifile,'tstart','long_name','Start processing day');
ncwriteatt(inifile,'tstart','units','day');

nccreate(inifile,'tend','Dimensions',{'one',1},'datatype','single');
ncwriteatt(inifile,'tend','long_name','End processing day');
ncwriteatt(inifile,'tend','units','day');

nccreate(inifile,'theta_s','Dimensions',{'one',1},'datatype','single');
ncwriteatt(inifile,'theta_s','long_name','S-coordinate surface control parameter');
ncwriteatt(inifile,'theta_s','units','nondimensional');

nccreate(inifile,'theta_b','Dimensions',{'one',1},'datatype','single');
ncwriteatt(inifile,'theta_b','long_name','S-coordinate bottom control parameter');
ncwriteatt(inifile,'theta_b','units','nondimensional');

nccreate(inifile,'Tcline','Dimensions',{'one',1},'datatype','single');
ncwriteatt(inifile,'Tcline','long_name','S-coordinate surface/bottom layer width');
ncwriteatt(inifile,'Tcline','units','meter');

nccreate(inifile,'hc','Dimensions',{'one',1},'datatype','single');
ncwriteatt(inifile,'hc','long_name','S-coordinate parameter critical depth');
ncwriteatt(inifile,'hc','units','meter');

nccreate(inifile,'sc_r','Dimensions',{'s_rho',N},'datatype','single');
ncwriteatt(inifile,'sc_r','long_name','S-coordinate at RHO-points');
ncwriteatt(inifile,'sc_r','units','-');

nccreate(inifile,'Cs_r','Dimensions',{'s_rho',N},'datatype','single');
ncwriteatt(inifile,'Cs_r','long_name','S-coordinate stretching curves at RHO-points');
ncwriteatt(inifile,'Cs_r','units','-');

nccreate(inifile,'ocean_time','Dimensions',{'time',1},'datatype','single');
ncwriteatt(inifile,'ocean_time','long_name','time since initialization');
ncwriteatt(inifile,'ocean_time','units','second');

nccreate(inifile,'u','Dimensions',{'xi_u',L,'eta_u',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'u','long_name','u-flux component');
ncwriteatt(inifile,'u','units','meter second-1');
%
nccreate(inifile,'v','Dimensions',{'xi_v',Lp,'eta_v',M,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'v','long_name','v-flux component');
ncwriteatt(inifile,'v','units','meter second-1');
%
nccreate(inifile,'w','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_w',Np,'time',1},'datatype','single');
ncwriteatt(inifile,'w','long_name','w-flux component');
ncwriteatt(inifile,'w','units','meter second-1');
%
nccreate(inifile,'ubar','Dimensions',{'xi_u',L,'eta_u',Mp,'time',1},'datatype','single');
ncwriteatt(inifile,'ubar','long_name','vertically integrated u-flux component');
ncwriteatt(inifile,'ubar','units','meter second-1');
%
nccreate(inifile,'vbar','Dimensions',{'xi_v',Lp,'eta_v',M,'time',1},'datatype','single');
ncwriteatt(inifile,'vbar','long_name','vertically integrated v-flux component');
ncwriteatt(inifile,'vbar','units','meter second-1');
%
nccreate(inifile,'zeta','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'time',1},'datatype','single');
ncwriteatt(inifile,'zeta','long_name','free surface');
ncwriteatt(inifile,'zeta','units','meter');
%
nccreate(inifile,'temp','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'temp','long_name','potential temperature');
ncwriteatt(inifile,'temp','units','Celcius');
%
nccreate(inifile,'salt','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'salt','long_name','Salinity');
ncwriteatt(inifile,'salt','units','PSU');
%
%  Write global attributes
%
 ncwriteatt(inifile,'/','Title',['R2R initial file for' gridfile]);
 ncwriteatt(inifile,'/','Date',date);
%
%
% Write some variables
%
ncwrite(inifile,'tstart',1.0);
ncwrite(inifile,'tend',1.0);
ncwrite(inifile,'theta_s',chdscd.theta_s);
ncwrite(inifile,'theta_b',chdscd.theta_b);
ncwrite(inifile,'Tcline',chdscd.hc);
ncwrite(inifile,'hc',chdscd.hc);
ncwrite(inifile,'ocean_time',0);

[sc_r,Cs_r] = sigma_stretch(chdscd.theta_s,chdscd.theta_b,N,'r',3);
disp('WARNING : writting sigma stretch sc and Cs')
disp('WARNING : check Sigma coord type in sigma_stretch.m')
ncwrite(inifile,'sc_r',sc_r);
ncwrite(inifile,'Cs_r',Cs_r);

disp('Init : w set to 0')
ncwrite(inifile,'w',zeros(Lp,Mp,Np));


%{bgc_tracer_list ;
%
%    BGC_INI.file_bgc{1} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_n_annual_landfilled.nc' ;   % '/data/project1/data/WOA18/Nitrate/woa18_all_n*nc';
%    BGC_INI.file_bgc{2} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_o_seasonal_landfilled.nc' ;  %'/data/project1/data/WOA18/Oxygen/woa18_all_o*nc';
%    BGC_INI.file_bgc{3} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_p_annual_landfilled.nc' ; %'/data/project1/data/WOA18/Phosphate/woa18_all_p*nc';
%    BGC_INI.file_bgc{4} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/woa13_all_i_seasonal_landfilled.nc'; %'/data/project1/data/WOA18/Silicate/woa18_all_i*nc';
%    BGC_INI.file_bgc{5} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/ETHZ_L1_46.1k_Anth_Mon.pop.h.1958-2006_Fe_1981-2000_Feymonavg_rename_add_time_dimvar.nc';
%    BGC_INI.file_bgc{6} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/n2ofromnn.nc' ;
%    BGC_INI.file_bgc{7} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/GLODAPv2.2016b.TCO2.nc' ;
%    BGC_INI.file_bgc{8} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/GLODAPv2.2016b.TAlk.nc' ;
%    BGC_INI.file_bgc{9} = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/SeaWiFS_CHL_MO_climatology_9km_landfill.nc';
%    BGC_INI.file_bgc{10}= '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/n2o_satn2o_sat_mon_mmolpm3_fromWOA13_landfill.nc' ;
%    BGC_INI.file_bgc{11}= '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/dpco2_climatology_T09_rev_landfill.nc';
%    BGC_INI.file_bgc{12}= '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/pacmed_0p25_ini2000_5daysAVGsoda.nc';
%
%    % file to read for forcing if type 'i'; value to init  if 'v' or scale factor if 'b'
%    BGC_INI.bgc_frcini = {BGC_INI.file_bgc{1},BGC_INI.file_bgc{3},BGC_INI.file_bgc{4},BGC_INI.file_bgc{5},BGC_INI.file_bgc{2},BGC_INI.file_bgc{9} ...
%                  '3.375','0.675','1.35e-05','0.0675', ...
%                  '0.2025','0.0675','1.35e-06','0.0675',...
%                  '0.0375','0.0075','7.5e-07', ...
%                  BGC_INI.file_bgc{7}, BGC_INI.file_bgc{8},'1.35', ...
%                  '1','0.8','0.1','0.003','0.0001',...
%                  BGC_INI.file_bgc{6}, ' ', BGC_INI.file_bgc{10}, '0.0001',BGC_INI.file_bgc{11},BGC_INI.file_bgc{12},...
%                  ' ',' '};
%
%BGC_INI.bgc_tracer = {'NO3','PO4','SiO3','Fe','O2','CHLA', ...
%                  'SPC','SPCHL','SPFE','SPCACO3', ...
%                  'DIATC','DIATCHL','DIATFE','DIATSI',...
%                  'DIAZC','DIAZCHL','DIAZFE', ...
%                  'DIC_glodap','Alk_glodap','ZOOC', ...
%                  'DON','DONR','DOP','DOPR','DOFE'...
%                  'N2O', 'N2O_SIDEN', 'N2O_ATM', 'N2O_NEV','pCO2','basindx', ...
%                  'DIC','Alk'};
%
%    BGC_INI.bgc_frctype = {'i','i','i','i','i','i' ...
%                  's_6','s_6','s_6','s_6', ...
%                  's_6','s_6','s_6','s_6',...
%                  's_6','s_6','s_6', ...
%                  'i','i','s_6', ...
%                  'v','v','v','v','v'...
%                  'i', 'b_26', 'i', 'v', 'i', 'i', ...
%                  'b_30','b_30'};
%
%    for trc=length(BGC_INI.bgc_frctype)+1:length(bgctracers_list.name)
%        BGC_INI.bgc_frctype{trc} = 'i';
%    end
%
%for trc=1:length(BGC_INI.bgc_tracer)
%
%    ind_trc = find(strcmp(BGC_INI.bgc_tracer{trc},bgctracers_list.name)) ;
%
%    if isempty(ind_trc)
%       disp(['ERROR : tracer '  BGC_INI.bgc_tracer{trc} ' not found in tracer list']) ; error ;
%    elseif ( strcmp(BGC_INI.bgc_tracer{trc},'pCO2')==1 || strcmp(BGC_INI.bgc_tracer{trc},'basindx')==1 )
%        nccreate(inifile,BGC_INI.bgc_tracer{trc},'Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'time',1},'datatype','single');
%        ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'long_name',bgctracers_list.longname{ind_trc});
%        ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'units',bgctracers_list.units{ind_trc});
%        type = BGC_INI.bgc_frctype{trc} ;
%        if strcmp(type(1),'i')==1
%           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',BGC_INI.bgc_frcini{trc});
%        elseif strcmp(type(1),'v')==1
%           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
%           % Write some variables
%%            var = ones(1,N,Mp,Lp)*str2num(BGC_INI.bgc_frcini{trc}) ;
%           var = ones(Lp,Mp,1)*str2num(BGC_INI.bgc_frcini{trc}) ;
%           ncwrite(inifile,BGC_INI.bgc_tracer{trc},var) ;
%        elseif strcmp(type(1),'b')==1
%           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
%        elseif strcmp(type(1),'s')==1
%           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['scaled form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
%        end
%    else
%        nccreate(inifile,BGC_INI.bgc_tracer{trc},'Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
%        ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'long_name',bgctracers_list.longname{ind_trc});
%        ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'units',bgctracers_list.units{ind_trc});
%        type = BGC_INI.bgc_frctype{trc} ;
%        if strcmp(type(1),'i')==1
%           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',BGC_INI.bgc_frcini{trc});
%        elseif strcmp(type(1),'v')==1
%           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
%           % Write some variables
%%            var = ones(1,N,Mp,Lp)*str2num(BGC_INI.bgc_frcini{trc}) ;
%           var = ones(Lp,Mp,N,1)*str2num(BGC_INI.bgc_frcini{trc}) ;
%           ncwrite(inifile,BGC_INI.bgc_tracer{trc},var) ;
%        elseif strcmp(type(1),'b')==1
%           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
%        end
%    end
%
%end
%}

nccreate(inifile,'PO4','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'PO4','long_name','Phosphate');
ncwriteatt(inifile,'PO4','units','mMol P m-3');

nccreate(inifile,'NO3','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'NO3','long_name','Nitrate');
ncwriteatt(inifile,'NO3','units','mMol N m-3');

nccreate(inifile,'SiO3','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'SiO3','long_name','Silicate');
ncwriteatt(inifile,'SiO3','units','mMol Si m-3');

nccreate(inifile,'NH4','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'NH4','long_name','Ammonium');
ncwriteatt(inifile,'NH4','units','mMol N m-3');

nccreate(inifile,'Fe','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'Fe','long_name','Iron');
ncwriteatt(inifile,'Fe','units','mMol Fe m-3');

nccreate(inifile,'O2','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'O2','long_name','Oxygen');
ncwriteatt(inifile,'O2','units','mMol O2 m-3');

nccreate(inifile,'DIC','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DIC','long_name','Dissolved inorganic carbon');
ncwriteatt(inifile,'DIC','units','mMol C m-3');

nccreate(inifile,'Alk','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'Alk','long_name','Alkalinity');
ncwriteatt(inifile,'Alk','units','mMol m-3');

nccreate(inifile,'DOC','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DOC','long_name','Dissolved organic carbon');
ncwriteatt(inifile,'DOC','units','mMol C m-3');

nccreate(inifile,'DON','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DON','long_name','Dissolved organic nitrogen');
ncwriteatt(inifile,'DON','units','mMol N m-3');

nccreate(inifile,'DOFE','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DOFE','long_name','Dissolved organic iron');
ncwriteatt(inifile,'DOFE','units','mMol Fe m-3');

nccreate(inifile,'DOP','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DOP','long_name','Dissolved organic phosphorus');
ncwriteatt(inifile,'DOP','units','mMol P m-3');

nccreate(inifile,'DOPR','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DOPR','long_name','Refractory dissolved organic phosphorus');
ncwriteatt(inifile,'DOPR','units','mMol P m-3');

nccreate(inifile,'DONR','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DONR','long_name','Refractory dissolved organic nitrogen');
ncwriteatt(inifile,'DONR','units','mMol N m-3');

nccreate(inifile,'ZOOC','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'ZOOC','long_name','Zooplankton');
ncwriteatt(inifile,'ZOOC','units','mMol C m-3');

nccreate(inifile,'SPC','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'SPC','long_name','Small phytoplankton carbon');
ncwriteatt(inifile,'SPC','units','mMol C m-3');

nccreate(inifile,'SPCHL','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'SPCHL','long_name','Small phytoplankton chlorophyll');
ncwriteatt(inifile,'SPCHL','units','mMol Chl-a m-3');

nccreate(inifile,'SPFE','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'SPFE','long_name','Small phytoplankton iron');
ncwriteatt(inifile,'SPFE','units','mMol Fe m-3');

nccreate(inifile,'SPCACO3','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'SPCACO3','long_name','Small phytoplankton CaCO3');
ncwriteatt(inifile,'SPCACO3','units','mMol CaCO3 m-3');

nccreate(inifile,'DIATC','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DIATC','long_name','Diatom carbon');
ncwriteatt(inifile,'DIATC','units','mMol C m-3');

nccreate(inifile,'DIATCHL','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DIATCHL','long_name','Diatom chlorophyll');
ncwriteatt(inifile,'DIATCHL','units','mMol Chl-a m-3');

nccreate(inifile,'DIATFE','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DIATFE','long_name','Diatom Iron');
ncwriteatt(inifile,'DIATFE','units','mMol Fe m-3');

nccreate(inifile,'DIATSI','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DIATSI','long_name','Diatom silicon');
ncwriteatt(inifile,'DIATSI','units','mMol Si m-3');

nccreate(inifile,'DIAZC','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DIAZC','long_name','Diazotroph carbon');
ncwriteatt(inifile,'DIAZC','units','mMol C m-3');

nccreate(inifile,'DIAZCHL','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DIAZCHL','long_name','Diazotroph chlorophyll');
ncwriteatt(inifile,'DIAZCHL','units','mMol Chl-a m-3');

nccreate(inifile,'DIAZFE','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'DIAZFE','long_name','Diazotroph iron');
ncwriteatt(inifile,'DIAZFE','units','mMol Fe m-3');

nccreate(inifile,'NO2','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'NO2','long_name','Nitrite');
ncwriteatt(inifile,'NO2','units','mMol N m-3');

nccreate(inifile,'N2','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'N2','long_name','Dinitrogen');
ncwriteatt(inifile,'N2','units','mMol N2 m-3');

nccreate(inifile,'N2O','Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
ncwriteatt(inifile,'N2O','long_name','Nitrous oxide');
ncwriteatt(inifile,'N2O','units','mMol N2O m-3');


return


