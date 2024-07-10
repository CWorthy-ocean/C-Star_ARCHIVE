%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %vec_starting from 1995
% %               1995 1996 1997 1998 1999
% time_offset = [0 365  366  365  365  365 ...
%                  366  365  365  365  366  365  365  365  366  365 ...
%                  365  365  366  365  365  365  366  365  365  365  366];
% %vec_starting from 2000
 time_offset = [0 366 365 365 365 366 365 365 365 366 365 ...
                  365 365 366 365 365 365 366 365 365 365 366];

for yr=2012:1:2012

    ind= yr - 1999;

%%%%%%%%%%%%%%%%%%%%% USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
%

%  frc_dir = '/glade/scratch/bachman/ROMS_tools/setup_s2r_phys+bgc/1.Make_grid/';
%  frcname = [frc_dir, 'Wales0_frc_bgc.nc'];

%
%%%%%%%%%%%%%%%%%%% END USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%

disp(['Working on ' frcname])

%%% BGC
if 1
ncid = netcdf.open(frcname,'NC_WRITE');
varid = netcdf.inqVarID(ncid,'pco2_time');
netcdf.reDef(ncid)
netcdf.delAtt(ncid,varid,'cycle_length')
netcdf.close(ncid)
time = ncread( frcname , 'pco2_time' ) ;
time = time + sum(time_offset(1:ind)) - 15.21875;
ncwrite(frcname , 'pco2_time' , time ) ;

ncid = netcdf.open(frcname,'NC_WRITE');
varid = netcdf.inqVarID(ncid,'iron_time');
netcdf.reDef(ncid)
netcdf.delAtt(ncid,varid,'cycle_length')
netcdf.close(ncid)
time = ncread( frcname , 'iron_time' ) ;
time = time + sum(time_offset(1:ind))  - 15.21875;
ncwrite(frcname , 'iron_time' , time ) ;

ncid = netcdf.open(frcname,'NC_WRITE');
varid = netcdf.inqVarID(ncid,'dust_time');
netcdf.reDef(ncid)
netcdf.delAtt(ncid,varid,'cycle_length')
netcdf.close(ncid)
time = ncread( frcname , 'dust_time' ) ;
time = time + sum(time_offset(1:ind))  - 15.21875;
ncwrite(frcname , 'dust_time' , time ) ;

ncid = netcdf.open(frcname,'NC_WRITE');
varid = netcdf.inqVarID(ncid,'NOx_time');
netcdf.reDef(ncid)
netcdf.delAtt(ncid,varid,'cycle_length')
netcdf.close(ncid)
time = ncread( frcname , 'NOx_time' ) ;
time = time + sum(time_offset(1:ind))  - 15.21875;
ncwrite(frcname , 'NOx_time' , time ) ;

ncid = netcdf.open(frcname,'NC_WRITE');
varid = netcdf.inqVarID(ncid,'NHy_time');
netcdf.reDef(ncid)
netcdf.delAtt(ncid,varid,'cycle_length')
netcdf.close(ncid)
time = ncread( frcname , 'NHy_time' ) ;
time = time + sum(time_offset(1:ind))  - 15.21875;
ncwrite(frcname , 'NHy_time' , time ) ;
end


end




