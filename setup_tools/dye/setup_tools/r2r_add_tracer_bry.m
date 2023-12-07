
bryname = '/glade/scratch/bachman/UCLA-ROMS/Work/Iceland4_dye/INPUT/ens8/Iceland4_bry.nc';

ncid = netcdf.open(bryname,'NC_WRITE');

s_rho_id = netcdf.inqDimID(ncid,'s_rho')
xi_rho_id = netcdf.inqDimID(ncid,'xi_rho')
eta_rho_id = netcdf.inqDimID(ncid,'eta_rho')
%bry_time_id = netcdf.inqDimID(ncid,'bry_time')
bry_time = ncread(bryname, 'bry_time');

[s_rho, s_rho_len] = netcdf.inqDim(ncid,s_rho_id)
[xi_rho, xi_rho_len] = netcdf.inqDim(ncid,xi_rho_id)
[eta_rho, eta_rho_len] = netcdf.inqDim(ncid,eta_rho_id)
%[bry_time, bry_time_len] = netcdf.inqDim(ncid,bry_time_id)
bry_time_len = numel(bry_time);

nccreate(bryname,'dye_south','Dimensions',{'xi_rho',xi_rho_len,'s_rho',s_rho_len,'time',bry_time_len},'datatype','single');
ncwriteatt(bryname,'dye_south','long_name','dye concentration along southern boundary');
ncwriteatt(bryname,'dye_south','units','nondimensional');

nccreate(bryname,'dye_east','Dimensions',{'eta_rho',eta_rho_len,'s_rho',s_rho_len,'time',bry_time_len},'datatype','single');
ncwriteatt(bryname,'dye_east','long_name','dye concentration along eastern boundary');
ncwriteatt(bryname,'dye_east','units','nondimensional');

nccreate(bryname,'dye_north','Dimensions',{'xi_rho',xi_rho_len,'s_rho',s_rho_len,'time',bry_time_len},'datatype','single');
ncwriteatt(bryname,'dye_north','long_name','dye concentration along northern boundary');
ncwriteatt(bryname,'dye_north','units','nondimensional');

nccreate(bryname,'dye_west','Dimensions',{'eta_rho',eta_rho_len,'s_rho',s_rho_len,'time',bry_time_len},'datatype','single');
ncwriteatt(bryname,'dye_west','long_name','dye concentration along western boundary');
ncwriteatt(bryname,'dye_west','units','nondimensional');

dye_south = zeros(xi_rho_len, s_rho_len, bry_time_len);
dye_east = zeros(eta_rho_len, s_rho_len, bry_time_len);
dye_north = zeros(xi_rho_len, s_rho_len, bry_time_len);
dye_west = zeros(eta_rho_len, s_rho_len, bry_time_len);

ncwrite(bryname, 'dye_south', dye_south);
ncwrite(bryname, 'dye_east', dye_east);
ncwrite(bryname, 'dye_north', dye_north);
ncwrite(bryname, 'dye_west', dye_west);

netcdf.close(ncid)
