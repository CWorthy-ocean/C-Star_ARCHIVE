function  create_grid(nx,ny,grdname,title)
%
% This is part of Easy Grid
%  (c) 2008, Jeroen Molemaker, UCLA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  Create dimensions
%

%
%  Create variables and attributes
%

%nw{'spherical'} = ncchar('one');
%nw{'spherical'}.long_name = 'Grid type logical switch';
%nw{'spherical'}.option_T = 'spherical';
nccreate(grdname,'spherical','Datatype','char','Dimensions', {'one',1});
ncwriteatt(grdname,'spherical','Long_name','Grid type logical switch');
ncwriteatt(grdname,'spherical','option_T','spherical');

%nw{'angle'} = ncdouble('eta_rho', 'xi_rho');
%nw{'angle'}.long_name = 'angle between xi axis and east';
%nw{'angle'}.units = 'radians';
nccreate(grdname,'angle','Dimensions', {'xi_rho',nx,'eta_rho',ny});
ncwriteatt(grdname,'angle','Long_name','Angle between xi axis and east');
ncwriteatt(grdname,'angle','units','radians');

%nw{'h'} = ncdouble('eta_rho', 'xi_rho');
%nw{'h'}.long_name = 'Final bathymetry at RHO-points';
%nw{'h'}.units = 'meter';
nccreate(grdname,'h','Dimensions', {'xi_rho',nx,'eta_rho',ny});
ncwriteatt(grdname,'h','Long_name','Final bathymetry at rho-points');
ncwriteatt(grdname,'h','units','meter');

%nw{'hraw'} = ncdouble('eta_rho', 'xi_rho');
%nw{'hraw'}.long_name = 'Working bathymetry at RHO-points';
%nw{'hraw'}.units = 'meter';
nccreate(grdname,'hraw','Dimensions', {'xi_rho',nx,'eta_rho',ny});
ncwriteatt(grdname,'hraw','Long_name','Working bathymetry at rho-points');
ncwriteatt(grdname,'hraw','units','meter');

%nw{'f'} = ncdouble('eta_rho', 'xi_rho');
%nw{'f'}.long_name = 'Coriolis parameter at RHO-points';
%nw{'f'}.units = 'second-1';
nccreate(grdname,'f','Dimensions', {'xi_rho',nx,'eta_rho',ny});
ncwriteatt(grdname,'f','Long_name','Coriolis parameter at rho-points');
ncwriteatt(grdname,'f','units','second-1');

%nw{'pm'} = ncdouble('eta_rho', 'xi_rho');
%nw{'pm'}.long_name = 'curvilinear coordinate metric in XI';
%nw{'pm'}.units = 'meter-1';
nccreate(grdname,'pm','Dimensions', {'xi_rho',nx,'eta_rho',ny});
ncwriteatt(grdname,'pm','Long_name','curvilinear coordinate metric in xi-direction');
ncwriteatt(grdname,'pm','units','meter-1');

%nw{'pn'} = ncdouble('eta_rho', 'xi_rho');
%nw{'pn'}.long_name = 'curvilinear coordinate metric in ETA';
%nw{'pn'}.units = 'meter-1';
nccreate(grdname,'pn','Dimensions', {'xi_rho',nx,'eta_rho',ny});
ncwriteatt(grdname,'pn','Long_name','curvilinear coordinate metric in eta-direction');
ncwriteatt(grdname,'pn','units','meter-1');

%nw{'lon_rho'} = ncdouble('eta_rho', 'xi_rho');
%nw{'lon_rho'}.long_name = 'longitude of RHO-points';
%nw{'lon_rho'}.units = 'degree_east';
nccreate(grdname,'lon_rho','Dimensions', {'xi_rho',nx,'eta_rho',ny});
ncwriteatt(grdname,'lon_rho','Long_name','longitude of rho-points');
ncwriteatt(grdname,'lon_rho','units','degree East');

%nw{'lat_rho'} = ncdouble('eta_rho', 'xi_rho');
%nw{'lat_rho'}.long_name = 'latitude of RHO-points';
%nw{'lat_rho'}.units = 'degree_north';
nccreate(grdname,'lat_rho','Dimensions', {'xi_rho',nx,'eta_rho',ny});
ncwriteatt(grdname,'lat_rho','Long_name','latitude of rho-points');
ncwriteatt(grdname,'lat_rho','units','degree North');

%nw{'mask_rho'} = ncdouble('eta_rho', 'xi_rho');
%nw{'mask_rho'}.long_name = 'mask on RHO-points';
%nw{'mask_rho'}.option_0 = 'land';
%nw{'mask_rho'}.option_1 = 'water';
nccreate(grdname,'mask_rho','Dimensions', {'xi_rho',nx,'eta_rho',ny});
ncwriteatt(grdname,'mask_rho','Long_name','mask at rho-points');
ncwriteatt(grdname,'mask_rho','units','land/water (0/1)');

%nw{'tra_lon'} = ncdouble('one');
%nw{'tra_lon'}.long_name = 'Easy grid: Longitudinal translation of base grid';
%nw{'tra_lon'}.units = 'degree East';
nccreate(grdname,'tra_lon','Dimensions', {'one',1});
ncwriteatt(grdname,'tra_lon','Long_name','Easy grid: Longitudinal translation of base grid');
ncwriteatt(grdname,'tra_lon','units','degree East');

%nw{'tra_lat'} = ncdouble('one');
%nw{'tra_lat'}.long_name = 'Easy grid: Latitudinal translation of base grid';
%nw{'tra_lat'}.units = 'degree North';
nccreate(grdname,'tra_lat','Dimensions', {'one',1});
ncwriteatt(grdname,'tra_lat','Long_name','Easy grid: Latitudinal translation of base grid');
ncwriteatt(grdname,'tra_lat','units','degree North');

%nw{'rotate'} = ncdouble('one');
%nw{'rotate'}.long_name = 'Easy grid: Rotation of base grid';
%nw{'rotate'}.units = 'degree';
nccreate(grdname,'rotate','Dimensions', {'one',1});
ncwriteatt(grdname,'rotate','Long_name','Easy grid: Rotation of base grid');
ncwriteatt(grdname,'rotate','units','degree');

%nw{'xy_flip'} = ncint('one');
%nw{'xy_flip'}.long_name = 'Easy grid: XY flip of base grid';
%nw{'xy_flip'}.units = 'True/False (0/1)';
nccreate(grdname,'xy_flip','Dimensions', {'one',1});
ncwriteatt(grdname,'xy_flip','Long_name','Easy grid: XY flip of base grid');
ncwriteatt(grdname,'xy_flip','units','True/False (0/1)');

%result = endef(nw);

%nw.title = title;
%nw.date = date;
%nw.type = 'ROMS grid produced by Easy Grid';
%nw.VertCoordType = 'NEW';
%result = close(nw);

ncwriteatt(grdname,'/','Title',title);
ncwriteatt(grdname,'/','Date',date);
ncwriteatt(grdname,'/','Type','ROMS grid produced by Easy Grid');

