function r2r_create_ini(inifile,gridfile,N,chdscd,makebgc,BGC_INI,bgctracers_list)
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

%{
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
ncwrite(inifile,'ocean_time',1.0*24*3600);

[sc_r,Cs_r] = sigma_stretch(chdscd.theta_s,chdscd.theta_b,N,'r',3);
disp('WARNING : writting sigma stretch sc and Cs')
disp('WARNING : check Sigma coord type in sigma_stretch.m')
ncwrite(inifile,'sc_r',sc_r);
ncwrite(inifile,'Cs_r',Cs_r);

disp('Init : w set to 0')
ncwrite(inifile,'w',zeros(Lp,Mp,Np));
%}

if makebgc==1

    for trc=1:length(BGC_INI.bgc_tracer)

    ind_trc = find(strcmp(BGC_INI.bgc_tracer{trc},bgctracers_list.name)) ;

    if isempty(ind_trc)
       disp(['ERROR : tracer '  BGC_INI.bgc_tracer{trc} ' not found in tracer list']) ; error ;
    elseif ( strcmp(BGC_INI.bgc_tracer{trc},'pCO2')==1 || strcmp(BGC_INI.bgc_tracer{trc},'basindx')==1 )
        nccreate(inifile,BGC_INI.bgc_tracer{trc},'Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'time',1},'datatype','single');
        ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
%            var = ones(1,N,Mp,Lp)*str2num(BGC_INI.bgc_frcini{trc}) ;
           var = ones(Lp,Mp,1)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(inifile,BGC_INI.bgc_tracer{trc},var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        elseif strcmp(type(1),'s')==1
           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['scaled form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end

        init_zero = zeros(Lp,Mp,1);
        ncwrite(inifile, BGC_INI.bgc_tracer{trc}, init_zero);
    else
        nccreate(inifile,BGC_INI.bgc_tracer{trc},'Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
        ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'long_name',bgctracers_list.longname{ind_trc});
        ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'units',bgctracers_list.units{ind_trc});
        type = BGC_INI.bgc_frctype{trc} ;
        if strcmp(type(1),'i')==1
           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',BGC_INI.bgc_frcini{trc});
        elseif strcmp(type(1),'v')==1
           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['homogeneous small concentration : ' BGC_INI.bgc_frcini{trc}]);
           % Write some variables
%            var = ones(1,N,Mp,Lp)*str2num(BGC_INI.bgc_frcini{trc}) ;
           var = ones(Lp,Mp,N,1)*str2num(BGC_INI.bgc_frcini{trc}) ;
           ncwrite(inifile,BGC_INI.bgc_tracer{trc},var) ;
        elseif strcmp(type(1),'b')==1
           ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'source',['build form : ' BGC_INI.bgc_tracer{str2num(type(3:end))}]);
        end

        init_zero = zeros(Lp,Mp,N,1);
        ncwrite(inifile, BGC_INI.bgc_tracer{trc}, init_zero);

    end

    end

end

return


