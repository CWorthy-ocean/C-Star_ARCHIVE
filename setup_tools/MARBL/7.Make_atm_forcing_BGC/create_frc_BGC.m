function create_frc_BGC(gridfile,frcfile,coarse_frc,vars)
%
%   Create ROMS bulk forcing file
%
%
if coarse_frc
  [nx,ny] = size(ncread(gridfile,'h_coarse'));
else
  [nx,ny] = size(ncread(gridfile,'h'));
end

if strcmp(vars,'pco2')

nccreate(frcfile,'pco2_time','Dimensions',{'time',12},'datatype','double');
ncwriteatt(frcfile,'pco2_time','long_name','time since 2000-1-1');
ncwriteatt(frcfile,'pco2_time','units','day');
ncwriteatt(frcfile,'pco2_time','cycle_length',365.25);

nccreate(frcfile,'pco2_air','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',12},'datatype','single');
ncwriteatt(frcfile,'pco2_air','long_name','atmospheric pco2');
ncwriteatt(frcfile,'pco2_air','units','ppmv');

elseif strcmp(vars,'dust')

nccreate(frcfile,'dust_time','Dimensions',{'time',12},'datatype','double');
ncwriteatt(frcfile,'dust_time','long_name','time since 2000-1-1');
ncwriteatt(frcfile,'dust_time','units','day');
ncwriteatt(frcfile,'dust_time','cycle_length',365.25);

nccreate(frcfile,'dust','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',12},'datatype','single');
ncwriteatt(frcfile,'dust','long_name','dust deposition');
ncwriteatt(frcfile,'dust','units','nmol/cm2/s');

elseif strcmp(vars,'iron')

nccreate(frcfile,'iron_time','Dimensions',{'time',12},'datatype','double');
ncwriteatt(frcfile,'iron_time','long_name','time since 2000-1-1');
ncwriteatt(frcfile,'iron_time','units','day');
ncwriteatt(frcfile,'iron_time','cycle_length',365.25);

nccreate(frcfile,'iron','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',12},'datatype','single');
ncwriteatt(frcfile,'iron','long_name','iron deposition');
ncwriteatt(frcfile,'iron','units','nmol/cm2/s');

elseif strcmp(vars,'NOx')

nccreate(frcfile,'NOx_time','Dimensions',{'time',12},'datatype','double');
ncwriteatt(frcfile,'NOx_time','long_name','days since 0001-01-01');
ncwriteatt(frcfile,'NOx_time','units','day');
ncwriteatt(frcfile,'NOx_time','cycle_length',365.25);

nccreate(frcfile,'NOx','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',12},'datatype','single');
ncwriteatt(frcfile,'NOx','long_name','NOx deposition');
ncwriteatt(frcfile,'NOx','units','g(N)/m2/s');

elseif strcmp(vars,'NHy')

nccreate(frcfile,'NHy_time','Dimensions',{'time',12},'datatype','double');
ncwriteatt(frcfile,'NHy_time','long_name','days since 0001-01-01');
ncwriteatt(frcfile,'NHy_time','units','day');
ncwriteatt(frcfile,'NHy_time','cycle_length',365.25);

nccreate(frcfile,'NHy','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',12},'datatype','single');
ncwriteatt(frcfile,'NHy','long_name','NHy deposition');
ncwriteatt(frcfile,'NHy','units','g(N)/m2/s');

else

disp('STOP : unable to create this variable')
STOP

end

%
%
%  Write global attributes
%
 ncwriteatt(frcfile,'/','Title','ROMS srf bgc forcings field');
 ncwriteatt(frcfile,'/','Date',date);
 ncwriteatt(frcfile,'/','gridfile',gridfile);
%
%
return


