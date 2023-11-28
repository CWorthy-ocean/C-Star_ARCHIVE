function create_bry(bryname,grdname,obcflag,param,BRYtime,BGC_INI,bgctracers_list);
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
T = length(BRYtime.source) ;

%
%  Create the boundary file
%

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
t_id    = netcdf.defDim(fid,'bry_time',length(BRYtime.source));
one_id  = netcdf.defDim(fid,'one',1);

netcdf.close(fid)

%
%  Create variables and attributes
nccreate(bryname,'theta_s','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'theta_s','long_name','S-coordinate surface control parameter');
ncwriteatt(bryname,'theta_s','units','nondimensional');

nccreate(bryname,'theta_b','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'theta_b','long_name','S-coordinate bottom control parameter');
ncwriteatt(bryname,'theta_b','units','nondimensional');

nccreate(bryname,'hc','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'hc','long_name','S-coordinate parameter critical depth');
ncwriteatt(bryname,'hc','units','meter');

nccreate(bryname,'Tcline','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'Tcline','long_name','S-coordinate surface/bottom layer width');
ncwriteatt(bryname,'Tcline','units','meter');

nccreate(bryname,'bry_time','Dimensions',{'bry_time',T},'datatype','single');
ncwriteatt(bryname,'bry_time','long_name','time for boundary data');
ncwriteatt(bryname,'bry_time','units','days');

nccreate(bryname,'tstart','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'tstart','long_name','start processing day');
ncwriteatt(bryname,'tstart','units','day');

nccreate(bryname,'tend','Dimensions',{'one',1},'datatype','single');
ncwriteatt(bryname,'tend','long_name','end processing day');
ncwriteatt(bryname,'tend','units','day');

nccreate(bryname,'sc_r','Dimensions',{'s_rho',N},'datatype','single');
ncwriteatt(bryname,'sc_r','long_name','S-coordinate at RHO-points');
ncwriteatt(bryname,'sc_r','units','-');

nccreate(bryname,'Cs_r','Dimensions',{'s_rho',N},'datatype','single');
ncwriteatt(bryname,'Cs_r','long_name','S-coordinate stretching curves at RHO-points');
ncwriteatt(bryname,'Cs_r','units','-');
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
ncwrite(bryname,'bry_time',[1:length(BRYtime.source)]);
ncwrite(bryname,'theta_s',theta_s);
ncwrite(bryname,'theta_b',theta_b);
ncwrite(bryname,'Tcline',hc);
ncwrite(bryname,'hc',hc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BGC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%a info = ncinfo(bryname);
% nvars = length(info.Variables);
%info.Variables.Name
% present = 0;
% for i=1:nvars
%   if strcmp(info.Variables(i).Name,'NO3_south')
%     present = 1;
%   end
% end

%  if ~present % can't add if already added

    for trc=1:length(BGC_INI.bgc_tracer)

    ind_trc = find(strcmp(BGC_INI.bgc_tracer{trc},bgctracers_list.name)) ;

    if isempty(ind_trc)

        disp(['ERROR : tracer '  BGC_INI.bgc_tracer{trc} ' not found in tracer list']); % error ;

    else   %%% 3D fields

        if obcflag(1)==1  %%   Southern boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_south'];
        nccreate(bryname,name_trc,'Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
 %       type = BGC_INI.bgc_frctype{trc} ;
 %       if strcmp(type(1),'i')==1
 %          ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
 %       elseif strcmp(type(1),'v')==1
 %          ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
 %          % Write some variables
 %          var = ones(Lp,N,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
 %          ncwrite(bryname,name_trc,var) ;
 %       elseif strcmp(type(1),'b')==1
 %          ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
 %       end
        init_zero = zeros(Lp,N,T);
        ncwrite(bryname,name_trc, init_zero);
        end

        if obcflag(2)==1  %%   Eastern boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_east'];
        nccreate(bryname,name_trc,'Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
 %       type = BGC_INI.bgc_frctype{trc} ;
 %       if strcmp(type(1),'i')==1
 %          ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
 %       elseif strcmp(type(1),'v')==1
 %          ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
 %          % Write some variables
 %          var = ones(Mp,N,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
 %          ncwrite(bryname,name_trc,var) ;
 %       elseif strcmp(type(1),'b')==1
 %          ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
 %       end
        init_zero = zeros(Mp,N,T);
        ncwrite(bryname,name_trc, init_zero);
        end

        if obcflag(3)==1  %%   Northern boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_north'];
        nccreate(bryname,name_trc,'Dimensions',{'xi_rho',Lp,'s_rho',N,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
 %       type = BGC_INI.bgc_frctype{trc} ;
 %       if strcmp(type(1),'i')==1
 %          ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
 %       elseif strcmp(type(1),'v')==1
 %          ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
 %          % Write some variables
 %          var = ones(Lp,N,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
 %          ncwrite(bryname,name_trc,var) ;
 %       elseif strcmp(type(1),'b')==1
 %          ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
 %       end
        init_zero = zeros(Lp,N,T);
        ncwrite(bryname,name_trc, init_zero);
        end

        if obcflag(4)==1  %%   Western boundary
        name_trc = [BGC_INI.bgc_tracer{trc} '_west'];
        nccreate(bryname,name_trc,'Dimensions',{'eta_rho',Mp,'s_rho',N,'bry_time',T},'datatype','single');
        ncwriteatt(bryname,name_trc,'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(bryname,name_trc,'units',bgctracers_list.units{ind_trc});
 %       type = BGC_INI.bgc_frctype{trc} ;
 %       if strcmp(type(1),'i')==1
 %          ncwriteatt(bryname,name_trc,'source',BGC_INI.bgc_frcini{trc});
 %       elseif strcmp(type(1),'v')==1
 %          ncwriteatt(bryname,name_trc,'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
 %          % Write some variables
 %          var = ones(Mp,N,T)*str2num(BGC_INI.bgc_frcini{trc}) ;
 %          ncwrite(bryname,name_trc,var) ;
 %       elseif strcmp(type(1),'b')==1
 %          ncwriteatt(bryname,name_trc,'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
 %       end
        init_zero = zeros(Mp,N,T);
        ncwrite(bryname,name_trc, init_zero);
        end

    end

    end

%  end  % if present

%end









