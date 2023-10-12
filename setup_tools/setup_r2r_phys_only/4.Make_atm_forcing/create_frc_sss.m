function create_frc_sss(gridfile,frcfile,coarse_frc)
%
%   Create ROMS bulk forcing file
% 
%
if coarse_frc
  [nx,ny] = size(ncread(gridfile,'h_coarse'));
else
  [nx,ny] = size(ncread(gridfile,'h'));
end

%
%  Create variables
%
nccreate(frcfile,'sss_time','Dimensions',{'time',12},'datatype','double');
ncwriteatt(frcfile,'sss_time','long_name','time since 2000-1-1');
ncwriteatt(frcfile,'sss_time','units','day');
ncwriteatt(frcfile,'sss_time','cycle_length',365.25);

nccreate(frcfile,'sss','Dimensions',{'xi_rho',nx,'eta_rho',ny,'time',12},'datatype','single');
ncwriteatt(frcfile,'sss','long_name','Sea surface Salinity');
ncwriteatt(frcfile,'sss','units','PSU');
ncwriteatt(frcfile,'sss','note','From WOA; 1955-2017');


%
%
%  Write global attributes
%
 ncwriteatt(frcfile,'/','Title','ROMS SSS field');
 ncwriteatt(frcfile,'/','Date',date);
 ncwriteatt(frcfile,'/','gridfile',gridfile);
%
%
return


