%  Generates the i,j locations of data extraction objects
%  for use as boundary forcing in a subsequent nested grid
%
%  Writes to netcdf file the i and j locations of the places 
%  from where we want to save data. Index locations are in [0,nx], [0,ny]
%  
%  It writes a data extraction object for each set of i,j points; each boundary
%  has 3 different objects, for rho,u, and v-points. The velocity objects
%  also include and angle to which the desired velocties will be rotated
%
%  If the child grid point is not in the parent domain it is given a value
%  of -1e5 
%
%  Note the mod statements for lonc and lonp. This is an attempt to deal
%  with parent and child longitudes that are possibly 360 apart
%  It will fail for grids that straddle the dateline of the zero meridian
%  In that case, subtract 180 degrees first
%

% -- START USER INPUT ----------
% Parent grid directory and file name
pdir    = '/paracas/nmolem/TWAKE/';
pname   = 'twake_grd.nc';
ename   = 'twake_edata.nc';

gname = 'SEB1';
lon = 121.4302;
lat = 22.7702;

gname = 'SEB2';
lon = 121.4428;
lat = 22.7978

% Output file name and info
info  = ['indices for ' gname ' in ' pname];


%lon =-122.507; 
%lat =  33.460;
%gname = 'CCE1'

%lon =-120.811; 
%lat =  34.305;
%gname = 'CCE2'

%lon =-122; 
%lat =  36.75;
%gname = 'MBARI_M1'

 period =  300;
 ang = 0;
 mooring_vars = 'zeta, temp, salt, u, v' ;

% -- END USER INPUT ------------

pname = [pdir pname];
ename = [pdir ename];

lonp = ncread(pname,'lon_rho');        
latp = ncread(pname,'lat_rho');
lonp = mod(lonp,360);

obj_name = gname;
obj_lon = lon; % + 360;
obj_lat = lat;
obj_ang = ang;
obj_msk =   1;
add_object(ename,obj_name,lonp,latp,period,obj_lon,obj_lat,obj_msk,obj_ang);
ncwriteatt(ename,obj_name,'output_vars',mooring_vars);

ncwriteatt(ename, '/', [gname '_info'],  info);           % info on parent and child grid




