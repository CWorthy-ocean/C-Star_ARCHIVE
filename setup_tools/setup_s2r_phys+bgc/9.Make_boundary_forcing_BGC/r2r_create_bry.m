function r2r_create_bry(bryname,grdname,obcflag,param,BRYtime,makebgc,BGC_INI,bgctracers_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   function l2r_create_bry(bryname,grdname,obcflag,...
%                          chdscd,cycle)
%
%   Input:
%
%   bryname      Netcdf climatology file name (character string)
%   grdname      Netcdf grid file name (character string)
%   obcflag      open boundary flag (1=open, [S E N W])
%   chdscd       S-coordinate parameters (object)
%   cycle        Length (days) for cycling the climatology (real)

%
%
% get S-coordinate parameters
%
theta_b = param.theta_b;
theta_s = param.theta_s;
hc      = param.hc;
N       = param.N;
%
%
%  Read the grid file and check the topography
%
h     = ncread(grdname, 'h')';
maskr = ncread(grdname, 'mask_rho')';
[Mp,Lp] = size(h);

hmin   = min(min(h(maskr==1)));

L  = Lp - 1;
M  = Mp - 1;
Np = N  + 1;
T = length(BRYtime.time) ;

%
%  Create the boundary file
%
%{

type    = 'BOUNDARY file';
history = 'ROMS';

nc_fileversion = netcdf.getConstant('NETCDF4');
fid = netcdf.create(bryname,nc_fileversion); % Here is the problem
netcdf.setFill(fid,'NOFILL');

%
%  Create dimensions
%
xi_uid  = netcdf.defDim(fid,'xi_u',L);
eta_vid = netcdf.defDim(fid,'eta_v',M);
xi_rid  = netcdf.defDim(fid,'xi_rho',Lp);
eta_rid = netcdf.defDim(fid,'eta_rho',Mp);
s_rid   = netcdf.defDim(fid,'s_rho',N);
s_wid   = netcdf.defDim(fid,'s_w',N+1);
t_id    = netcdf.defDim(fid,'bry_time',length(BRYtime.time));
one_id  = netcdf.defDim(fid,'one',1);

netcdf.close(fid)

%
%  Create variables and attributes
%
nccreate(bryname,'theta_s','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'theta_s','long_name','S-coordinate surface control parameter');
ncwriteatt(bryname,'theta_s','units','nondimensional');
%
nccreate(bryname,'theta_b','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'theta_b','long_name','S-coordinate bottom control parameter');
ncwriteatt(bryname,'theta_b','units','nondimensional');
%
nccreate(bryname,'hc','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'hc','long_name','S-coordinate parameter critical depth');
ncwriteatt(bryname,'hc','units','meter');
%
nccreate(bryname,'Tcline','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'Tcline','long_name','S-coordinate surface/bottom layer width');
ncwriteatt(bryname,'Tcline','units','meter');
%
nccreate(bryname,'bry_time','Dimensions',{'bry_time',T},'datatype','single');
ncwriteatt(bryname,'bry_time','long_name','time for boundary data');
ncwriteatt(bryname,'bry_time','units','days');
if BRYtime.cycle>0
 ncwriteatt(bryname,'bry_time','cycle_length',BRYtime.cycle);
end
%
nccreate(bryname,'tstart','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'tstart','long_name','start processing day');
ncwriteatt(bryname,'tstart','units','day');
%
nccreate(bryname,'tend','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'tend','long_name','end processing day');
ncwriteatt(bryname,'tend','units','day');
%
nccreate(bryname,'sc_r','Dimensions',{'s_rho',N},'datatype','single');
ncwriteatt(bryname,'sc_r','long_name','S-coordinate at RHO-points');
ncwriteatt(bryname,'sc_r','units','-');
%
nccreate(bryname,'Cs_r','Dimensions',{'s_rho',N},'datatype','single');
ncwriteatt(bryname,'Cs_r','long_name','S-coordinate stretching curves at RHO-points');
ncwriteatt(bryname,'Cs_r','units','-');
%
if obcflag(1)==1  %%   Southern boundary
%
%   nc{'temp_south'} = ncfloat('bry_time','s_rho','xi_rho') ;
%   nc{'temp_south'}.long_name = 'southern boundary potential temperature';
%   nc{'temp_south'}.units = 'Celsius';
  nccreate(bryname,'temp_south','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'temp_south','long_name','southern boundary potential temperature');
  ncwriteatt(bryname,'temp_south','units','Celsius');
%
%   nc{'salt_south'} = ncfloat('bry_time','s_rho','xi_rho') ;
%   nc{'salt_south'}.long_name = 'southern boundary salinity';
%   nc{'salt_south'}.units = 'PSU';
  nccreate(bryname,'salt_south','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'salt_south','long_name','southern boundary salinity');
  ncwriteatt(bryname,'salt_south','units','PSU');
%
%   nc{'u_south'} = ncfloat('bry_time','s_rho','xi_u') ;
%   nc{'u_south'}.long_name = 'southern boundary u-momentum component';
%   nc{'u_south'}.units = 'meter second-1';
  nccreate(bryname,'u_south','Dimensions',{'xi_u',L,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'u_south','long_name','southern boundary u-momentum component');
  ncwriteatt(bryname,'u_south','units','meter second-1');
%
%   nc{'v_south'} = ncfloat('bry_time','s_rho','xi_rho') ;
%   nc{'v_south'}.long_name = 'southern boundary v-momentum component';
%   nc{'v_south'}.units = 'meter second-1';
  nccreate(bryname,'v_south','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'v_south','long_name','southern boundary v-momentum component');
  ncwriteatt(bryname,'v_south','units','meter second-1');
%
%   nc{'ubar_south'} = ncfloat('bry_time','xi_u') ;
%   nc{'ubar_south'}.long_name = 'southern boundary vertically integrated u-momentum component';
%   nc{'ubar_south'}.units = 'meter second-1';
  nccreate(bryname,'ubar_south','Dimensions',{'xi_u',L,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'ubar_south','long_name','southern boundary vertically integrated u-momentum component');
  ncwriteatt(bryname,'ubar_south','units','meter second-1');
%
%   nc{'vbar_south'} = ncfloat('bry_time','xi_rho') ;
%   nc{'vbar_south'}.long_name = 'southern boundary vertically integrated v-momentum component';
%   nc{'vbar_south'}.units = 'meter second-1';
  nccreate(bryname,'vbar_south','Dimensions',{'xi_rho',Lp,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'vbar_south','long_name','southern boundary vertically integrated v-momentum component');
  ncwriteatt(bryname,'vbar_south','units','meter second-1');
%
%   nc{'zeta_south'} = ncfloat('bry_time','xi_rho') ;
%   nc{'zeta_south'}.long_name = 'southern boundary sea surface height';
%   nc{'zeta_south'}.units = 'meter';
  nccreate(bryname,'zeta_south','Dimensions',{'xi_rho',Lp,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'zeta_south','long_name','southern boundary sea surface height');
  ncwriteatt(bryname,'zeta_south','units','meter');
end
%
%
if obcflag(2)==1  %%   Eastern boundary
%
%   nc{'temp_east'} = ncfloat('bry_time','s_rho','eta_rho') ;
%   nc{'temp_east'}.long_name = 'eastern boundary potential temperature';
%   nc{'temp_east'}.units = 'Celsius';
  nccreate(bryname,'temp_east','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'temp_east','long_name','eastern boundary potential temperature');
  ncwriteatt(bryname,'temp_east','units','Celsius');
%
%   nc{'salt_east'} = ncfloat('bry_time','s_rho','eta_rho') ;
%   nc{'salt_east'}.long_name = 'eastern boundary salinity';
%   nc{'salt_east'}.units = 'PSU';
  nccreate(bryname,'salt_east','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'salt_east','long_name','eastern boundary salinity');
  ncwriteatt(bryname,'salt_east','units','PSU');
%
%   nc{'u_east'} = ncfloat('bry_time','s_rho','eta_rho') ;
%   nc{'u_east'}.long_name = 'eastern boundary u-momentum component';
%   nc{'u_east'}.units = 'meter second-1';
  nccreate(bryname,'u_east','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'u_east','long_name','eastern boundary u-momentum component');
  ncwriteatt(bryname,'u_east','units','meter second-1');
%
%   nc{'v_east'} = ncfloat('bry_time','s_rho','eta_v') ;
%   nc{'v_east'}.long_name = 'eastern boundary v-momentum component';
%   nc{'v_east'}.units = 'meter second-1';
  nccreate(bryname,'v_east','Dimensions',{'eta_v',M,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'v_east','long_name','eastern boundary v-momentum component');
  ncwriteatt(bryname,'v_east','units','meter second-1');
%
%   nc{'ubar_east'} = ncfloat('bry_time','eta_rho') ;
%   nc{'ubar_east'}.long_name = 'eastern boundary vertically integrated u-momentum component';
%   nc{'ubar_east'}.units = 'meter second-1';
  nccreate(bryname,'ubar_east','Dimensions',{'eta_rho',Mp,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'ubar_east','long_name','eastern boundary vertically integrated u-momentum component');
  ncwriteatt(bryname,'ubar_east','units','meter second-1');
%
%   nc{'vbar_east'} = ncfloat('bry_time','eta_v') ;
%   nc{'vbar_east'}.long_name = 'eastern boundary vertically integrated v-momentum component';
%   nc{'vbar_east'}.units = 'meter second-1';
  nccreate(bryname,'vbar_east','Dimensions',{'eta_v',M,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'vbar_east','long_name','eastern boundary vertically integrated v-momentum component');
  ncwriteatt(bryname,'vbar_east','units','meter second-1');
%
%   nc{'zeta_east'} = ncfloat('bry_time','eta_rho') ;
%   nc{'zeta_east'}.long_name = 'eastern boundary sea surface height';
%   nc{'zeta_east'}.units = 'meter';
  nccreate(bryname,'zeta_east','Dimensions',{'eta_rho',Mp,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'zeta_east','long_name','eastern boundary sea surface height');
  ncwriteatt(bryname,'zeta_east','units','meter');
end
%
if obcflag(3)==1  %%   Northern boundary
%
%   nc{'temp_north'} = ncfloat('bry_time','s_rho','xi_rho') ;
%   nc{'temp_north'}.long_name = 'northern boundary potential temperature';
%   nc{'temp_north'}.units = 'Celsius';
  nccreate(bryname,'temp_north','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'temp_north','long_name','northern boundary potential temperature');
  ncwriteatt(bryname,'temp_north','units','Celsius');
%
%   nc{'salt_north'} = ncfloat('bry_time','s_rho','xi_rho') ;
%   nc{'salt_north'}.long_name = 'northern boundary salinity';
%   nc{'salt_north'}.units = 'PSU';
  nccreate(bryname,'salt_north','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'salt_north','long_name','northern boundary salinity');
  ncwriteatt(bryname,'salt_north','units','PSU');
%
%   nc{'u_north'} = ncfloat('bry_time','s_rho','xi_u') ;
%   nc{'u_north'}.long_name = 'northern boundary u-momentum component';
%   nc{'u_north'}.units = 'meter second-1';
  nccreate(bryname,'u_north','Dimensions',{'xi_u',L,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'u_north','long_name','orthern boundary u-momentum component');
  ncwriteatt(bryname,'u_north','units','meter second-1');
%
%   nc{'v_north'} = ncfloat('bry_time','s_rho','xi_rho') ;
%   nc{'v_north'}.long_name = 'northern boundary v-momentum component';
%   nc{'v_north'}.units = 'meter second-1';
  nccreate(bryname,'v_north','Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'v_north','long_name','northern boundary v-momentum component');
  ncwriteatt(bryname,'v_north','units','meter second-1');
%
%   nc{'ubar_north'} = ncfloat('bry_time','xi_u') ;
%   nc{'ubar_north'}.long_name = 'northern boundary vertically integrated u-momentum component';
%   nc{'ubar_north'}.units = 'meter second-1';
  nccreate(bryname,'ubar_north','Dimensions',{'xi_u',L,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'ubar_north','long_name','northern boundary vertically integrated u-momentum component');
  ncwriteatt(bryname,'ubar_north','units','meter second-1');
%
%   nc{'vbar_north'} = ncfloat('bry_time','xi_rho') ;
%   nc{'vbar_north'}.long_name = 'northern boundary vertically integrated v-momentum component';
%   nc{'vbar_north'}.units = 'meter second-1';
  nccreate(bryname,'vbar_north','Dimensions',{'xi_rho',Lp,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'vbar_north','long_name','northern boundary vertically integrated v-momentum component');
  ncwriteatt(bryname,'vbar_north','units','meter second-1');
%
%   nc{'zeta_north'} = ncfloat('bry_time','xi_rho') ;
%   nc{'zeta_north'}.long_name = 'northern boundary sea surface height';
%   nc{'zeta_north'}.units = 'meter';
  nccreate(bryname,'zeta_north','Dimensions',{'xi_rho',Lp,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'zeta_north','long_name','northern boundary sea surface height');
  ncwriteatt(bryname,'zeta_north','units','meter');
end
%
if obcflag(4)==1  %%   Western boundary
%
%   nc{'temp_west'} = ncfloat('bry_time','s_rho','eta_rho') ;
%   nc{'temp_west'}.long_name = 'western boundary potential temperature';
%   nc{'temp_west'}.units = 'Celsius';
  nccreate(bryname,'temp_west','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'temp_west','long_name','western boundary potential temperature');
  ncwriteatt(bryname,'temp_west','units','Celsius');
%
%   nc{'salt_west'} = ncfloat('bry_time','s_rho','eta_rho') ;
%   nc{'salt_west'}.long_name = 'western boundary salinity';
%   nc{'salt_west'}.units = 'PSU';
  nccreate(bryname,'salt_west','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'salt_west','long_name','western boundary salinity');
  ncwriteatt(bryname,'salt_west','units','PSU');
%
%   nc{'u_west'} = ncfloat('bry_time','s_rho','eta_rho') ;
%   nc{'u_west'}.long_name = 'western boundary u-momentum component';
%   nc{'u_west'}.units = 'meter second-1';
  nccreate(bryname,'u_west','Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'u_west','long_name','western boundary u-momentum component');
  ncwriteatt(bryname,'u_west','units','meter second-1');
%
%   nc{'v_west'} = ncfloat('bry_time','s_rho','eta_v') ;
%   nc{'v_west'}.long_name = 'western boundary v-momentum component';
%   nc{'v_west'}.units = 'meter second-1';
  nccreate(bryname,'v_west','Dimensions',{'eta_v',M,'s_rho',N,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'v_west','long_name','western boundary v-momentum component');
  ncwriteatt(bryname,'v_west','units','meter second-1');
%
%   nc{'ubar_west'} = ncfloat('bry_time','eta_rho') ;
%   nc{'ubar_west'}.long_name = 'western boundary vertically integrated u-momentum component';
%   nc{'ubar_west'}.units = 'meter second-1';
  nccreate(bryname,'ubar_west','Dimensions',{'eta_rho',Mp,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'ubar_west','long_name','western boundary vertically integrated u-momentum component');
  ncwriteatt(bryname,'ubar_west','units','meter second-1');
%
%   nc{'vbar_west'} = ncfloat('bry_time','eta_v') ;
%   nc{'vbar_west'}.long_name = 'western boundary vertically integrated v-momentum component';
%   nc{'vbar_west'}.units = 'meter second-1';
  nccreate(bryname,'vbar_west','Dimensions',{'eta_v',M,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'vbar_west','long_name','western boundary vertically integrated v-momentum component');
  ncwriteatt(bryname,'vbar_west','units','meter second-1');
%
%   nc{'zeta_west'} = ncfloat('bry_time','eta_rho') ;
%   nc{'zeta_west'}.long_name = 'western boundary sea surface height';
%   nc{'zeta_west'}.units = 'meter';
  nccreate(bryname,'zeta_west','Dimensions',{'eta_rho',Mp,'bry_time',T},'datatype','single');
  ncwriteatt(bryname,'zeta_west','long_name','western boundary sea surface height');
  ncwriteatt(bryname,'zeta_west','units','meter');
end
%
%
% Create global attributes
%
 ncwriteatt(bryname,'/','Title',['Boundary file for' grdname]);
 ncwriteatt(bryname,'/','Date',date);
 ncwriteatt(bryname,'/','type',type);
 ncwriteatt(bryname,'/','history',history);

%
% Write variables
%
ncwrite(bryname,'tstart',BRYtime.tstart);
ncwrite(bryname,'tend',BRYtime.tend);
ncwrite(bryname,'bry_time',BRYtime.time);
ncwrite(bryname,'theta_s',theta_s);
ncwrite(bryname,'theta_b',theta_b);
ncwrite(bryname,'Tcline',hc);
ncwrite(bryname,'hc',hc);


[sc_r,Cs_r] = sigma_stretch(theta_s,theta_b,N,'r',3);
disp('WARNING : writting sigma stretch sc and Cs')
disp('WARNING : check Sigma coord type in sigma_stretch.m')
ncwrite(bryname,'sc_r',sc_r);
ncwrite(bryname,'Cs_r',Cs_r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BGC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if makebgc==1
%}

 info = ncinfo(bryname);
 nvars = length(info.Variables);
info.Variables.Name
 present = 0;
 for i=1:nvars
   if strcmp(info.Variables(i).Name,'NO3_south')
     present = 1;
   end
 end

  if ~present % can't add if already added

    for trc=1:length(BGC_INI.bgc_tracer)

    ind_trc = find(strcmp(BGC_INI.bgc_tracer{trc},bgctracers_list.name)) ;

    if isempty(ind_trc)

        disp(['ERROR : tracer '  BGC_INI.bgc_tracer{trc} ' not found in tracer list']) ; error ;

    elseif ( strcmp(BGC_INI.bgc_tracer{trc},'pCO2')==1 || strcmp(BGC_INI.bgc_tracer{trc},'basindx')==1 )

        if obcflag(1)==1  %%   Southern boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_south'];
        nccreate(bryname,name_trc,'Dimensions',{'xi_rho',Lp,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
           var = ones(Lp,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(bryname,name_trc,var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        elseif strcmp(type(1),'s')==1
           ncwriteatt(bryname,name_trc,'source',['scaled form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end
        end

        if obcflag(2)==1  %%   Eastern boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_east'];
        nccreate(bryname,name_trc,'Dimensions',{'eta_rho',Mp,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
           var = ones(Mp,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(bryname,name_trc,var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        elseif strcmp(type(1),'s')==1
           ncwriteatt(bryname,name_trc,'source',['scaled form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end
        end

        if obcflag(3)==1  %%   Northern boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_north'];
        nccreate(bryname,name_trc,'Dimensions',{'xi_rho',Lp,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
           var = ones(Lp,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(bryname,name_trc,var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        elseif strcmp(type(1),'s')==1
           ncwriteatt(bryname,name_trc,'source',['scaled form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end
        end

        if obcflag(4)==1  %%   Western boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_west'];
        nccreate(bryname,name_trc,'Dimensions',{'eta_rho',Mp,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
           var = ones(Mp,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(bryname,name_trc,var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        elseif strcmp(type(1),'s')==1
           ncwriteatt(bryname,name_trc,'source',['scaled form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end
        end

    else   %%% 3D fields

        if obcflag(1)==1  %%   Southern boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_south'];
        nccreate(bryname,name_trc,'Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
           var = ones(Lp,N,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(bryname,name_trc,var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end
        init_zero = zeros(Lp,N,T);
        ncwrite(bryname,name_trc, init_zero);
        end

        if obcflag(2)==1  %%   Eastern boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_east'];
        nccreate(bryname,name_trc,'Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
           var = ones(Mp,N,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(bryname,name_trc,var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end
        init_zero = zeros(Mp,N,T);
        ncwrite(bryname,name_trc, init_zero);
        end

        if obcflag(3)==1  %%   Northern boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_north'];
        nccreate(bryname,name_trc,'Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
           var = ones(Lp,N,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(bryname,name_trc,var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end
        init_zero = zeros(Lp,N,T);
        ncwrite(bryname,name_trc, init_zero);
        end

        if obcflag(4)==1  %%   Western boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_west'];
        nccreate(bryname,name_trc,'Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
           var = ones(Mp,N,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(bryname,name_trc,var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end
        init_zero = zeros(Mp,N,T);
        ncwrite(bryname,name_trc, init_zero);
        end

    end

    end

  end  % if present

end









