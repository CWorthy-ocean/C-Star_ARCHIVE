function create_ini(inifile,gridfile,N,chdscd,BGC_INI,bgctracers_list)

%  Read the grid file

h       = ncread(gridfile,'h')';
[Mp,Lp] = size(h);
L       = Lp - 1 ;
M       = Mp - 1 ;
Np      = N  + 1 ;

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

%  Write global attributes

ncwriteatt(inifile,'/','Title',['W2R initial file for' gridfile]);
ncwriteatt(inifile,'/','Date',date);

% Write some variables

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


for trc=1:length(BGC_INI.bgc_tracer)

  ind_trc = find(strcmp(BGC_INI.bgc_tracer{trc},bgctracers_list.name)) ;

  if isempty(ind_trc)
     disp(['ERROR : tracer '  BGC_INI.bgc_tracer{trc} ' not found in tracer list']) ; error ;
  else
     nccreate(inifile,BGC_INI.bgc_tracer{trc},'Dimensions',{'xi_rho',Lp,'eta_rho',Mp,'s_rho',N,'time',1},'datatype','single');
     ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'long_name',bgctracers_list.longname{ind_trc});
     ncwriteatt(inifile,BGC_INI.bgc_tracer{trc},'units',bgctracers_list.units{ind_trc});

     init_zero = zeros(Lp,Mp,N,1);
     ncwrite(inifile, BGC_INI.bgc_tracer{trc}, init_zero);

  end

end

return


