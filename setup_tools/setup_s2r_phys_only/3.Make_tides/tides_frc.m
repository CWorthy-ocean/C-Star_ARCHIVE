
% Re-factor netcdf files for tpxo data.

fname_in =  'DATA/grid_tpxo9.v2a.nc';

h = ncread(fname_in,'hz')';


fname_in =  'DATA/h_tpxo9.v2a.nc';

con = ncread(fname_in,'con');

lon = ncread(fname_in,'lon_z')';
lat = ncread(fname_in,'lat_z')';

hr = ncread(fname_in,'hRe');
hr = permute(hr,[2 1 3]);
hi = ncread(fname_in,'hIm');
hi = permute(hi,[2 1 3]);

fname_in =  'DATA/u_tpxo9.v2a.nc';

lon_u = ncread(fname_in,'lon_u')';
lat_u = ncread(fname_in,'lat_u')';
lon_v = ncread(fname_in,'lon_v')';
lat_v = ncread(fname_in,'lat_v')';

ur = ncread(fname_in,'URe');
ur = permute(ur,[2 1 3]);
ui = ncread(fname_in,'UIm');
ui = permute(ui,[2 1 3]);

vr = ncread(fname_in,'VRe');
vr = permute(vr,[2 1 3]);
vi = ncread(fname_in,'VIm');
vi = permute(vi,[2 1 3]);

fname_in =  'DATA/sal_tpxo9.v2a.nc';

sr = ncread(fname_in,'hRe');
sr = permute(sr,[2 1 3]);
si = ncread(fname_in,'hIm');
si = permute(si,[2 1 3]);

om = [1.405189e-04, 1.454441e-04, 1.378797e-04, 1.458423e-04, ...
      7.292117e-05, 6.759774e-05, 7.252295e-05, 6.495854e-05, ...
      0.026392e-04, 0.053234e-04, 2.810377e-04, 2.783984e-04, ...
      2.859630e-04, 1.352405e-04, 7.2722e-05 ];

fname_out = 'tpxo9.v2a.nc';

   nx = 2160 ;
   ny = 1081 ;
   nc =   15 ;

%  nccreate(fname_out,'constituents','dimensions',{'nc',nc},'datatype','char');
%  ncwriteatt(fname_out,'constituents','long_name','Tidal constituents');
%  ncwrite(fname_out,'constituents',con);

   nccreate(fname_out,'omega','dimensions',{'nc',nc},'datatype','double');
   ncwriteatt(fname_out,'omega','long_name','Tidal frequencies');
   ncwriteatt(fname_out,'omega','units','1/s');
   ncwrite(fname_out,'omega',om);

   %----- Lat/Lon----
   nccreate(fname_out,'lon_r','dimensions',{'nx',nx,'ny',ny},'datatype','single');
   ncwriteatt(fname_out,'lon_r','long_name','Longitude at rho points');
   ncwriteatt(fname_out,'lon_r','units','Degrees East');
   ncwrite(fname_out,'lon_r',lon);

   nccreate(fname_out,'lat_r','dimensions',{'nx',nx,'ny',ny},'datatype','single');
   ncwriteatt(fname_out,'lat_r','long_name','Latitude at rho points');
   ncwriteatt(fname_out,'lat_r','units','Degrees North');
   ncwrite(fname_out,'lat_r',lat);

   nccreate(fname_out,'lon_u','dimensions',{'nx',nx,'ny',ny},'datatype','single');
   ncwriteatt(fname_out,'lon_u','long_name','Longitude at u-points');
   ncwriteatt(fname_out,'lon_u','units','Degrees East');
   ncwrite(fname_out,'lon_u',lon_u);

   nccreate(fname_out,'lat_u','dimensions',{'nx',nx,'ny',ny},'datatype','single');
   ncwriteatt(fname_out,'lat_u','long_name','Latitude at u-points');
   ncwriteatt(fname_out,'lat_u','units','Degrees North');
   ncwrite(fname_out,'lat_u',lat_u);

   nccreate(fname_out,'lon_v','dimensions',{'nx',nx,'ny',ny},'datatype','single');
   ncwriteatt(fname_out,'lon_v','long_name','Longitude at v-points');
   ncwriteatt(fname_out,'lon_v','units','Degrees East');
   ncwrite(fname_out,'lon_v',lon_u);

   nccreate(fname_out,'lat_v','dimensions',{'nx',nx,'ny',ny},'datatype','single');
   ncwriteatt(fname_out,'lat_v','long_name','Latitude at v-points');
   ncwriteatt(fname_out,'lat_v','units','Degrees North');
   ncwrite(fname_out,'lat_v',lat_u);

   nccreate(fname_out,'depth','dimensions',{'nx',nx,'ny',ny},'datatype','single');
   ncwriteatt(fname_out,'depth','long_name','Bathymetry at rho-points');
   ncwriteatt(fname_out,'depth','units','meter');
   ncwrite(fname_out,'depth',h);

   %----- fields ----
   nccreate(fname_out,'h_Re','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
   ncwriteatt(fname_out,'h_Re','long_name','Tidal elevation, complex; real part');
   ncwriteatt(fname_out,'h_Re','units','meter');
   ncwrite(fname_out,'h_Re',hr);

   nccreate(fname_out,'h_Im','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
   ncwriteatt(fname_out,'h_Im','long_name','Tidal elevation, complex; imag part');
   ncwriteatt(fname_out,'h_Im','units','meter');
   ncwrite(fname_out,'h_Im',hi);

   nccreate(fname_out,'sal_Re','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
   ncwriteatt(fname_out,'sal_Re','long_name','SAL elevation, complex; real part');
   ncwriteatt(fname_out,'sal_Re','units','meter');
   ncwrite(fname_out,'sal_Re',sr);

   nccreate(fname_out,'sal_Im','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
   ncwriteatt(fname_out,'sal_Im','long_name','SAL elevation, complex; imag part');
   ncwriteatt(fname_out,'sal_Im','units','meter');
   ncwrite(fname_out,'sal_Im',si);

   nccreate(fname_out,'u_Re','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
   ncwriteatt(fname_out,'u_Re','long_name','Tidal transport WE, complex; real part');
   ncwriteatt(fname_out,'u_Re','units','meter^2/s');
   ncwrite(fname_out,'u_Re',ur);

   nccreate(fname_out,'u_Im','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
   ncwriteatt(fname_out,'u_Im','long_name','Tidal transport WE, complex; imag part');
   ncwriteatt(fname_out,'u_Im','units','meter^2/s');
   ncwrite(fname_out,'u_Im',ui);

   nccreate(fname_out,'v_Re','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
   ncwriteatt(fname_out,'v_Re','long_name','Tidal transport SN, complex; real part');
   ncwriteatt(fname_out,'v_Re','units','meter^2/s');
   ncwrite(fname_out,'v_Re',vr);

   nccreate(fname_out,'v_Im','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
   ncwriteatt(fname_out,'v_Im','long_name','Tidal transport SN, complex; imag part');
   ncwriteatt(fname_out,'v_Im','units','meter^2/s');
   ncwrite(fname_out,'v_Im',vi);


%  nccreate(fname_out,'h_re','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
%  ncwriteatt(fname_out,'h_re','long_name','Tidal elevation complex Relitude, real part');
%  ncwriteatt(fname_out,'h_re','units','meter');
%  ncwrite(fname_out,'h_re',ha);

%  nccreate(fname_out,'h_im','dimensions',{'nx',nx,'ny',ny, 'nc',nc},'datatype','single');
%  ncwriteatt(fname_out,'h_im','long_name','Tidal elevation complex Relitude, imag part');
%  ncwriteatt(fname_out,'h_im','units','meter');
%  ncwrite(fname_out,'h_im',ha);



   ncwriteatt(fname_out,'/','type','OTIS tidal elevation file and SAL');
   ncwriteatt(fname_out,'/','version','TPXO9.v2a 2020: deep=TPXO9.v1, shallow - averaged TPXO9-atlas-v2');
   ncwriteatt(fname_out,'/','Constituents',con);



%     data omega_d/ &
%       m2 1.405189e-04, 
%      	s2 1.454441e-04, 
%       n2 1.378797e-04, 
%       k2 1.458423e-04, 
%       k1 7.292117e-05,
%      	o1 6.759774e-05,
%       p1 7.252295e-05, 
%       q1 6.495854e-05,
%      	mm 0.026392e-04,
%      	mf 0.053234e-04,
%       m4  2.810377e-04,
%       mn4 2.783984e-04, 
%      	ms4 2.859630e-04, 
%       2n2 1.352405e-04,
%      	s1  7.2722e-05/
%      	mu2 1.355937e-04,
%      	nu2 1.382329e-04,
%      	l2  1.431581e-04,
%       t2  1.452450e-04,
%      	j1 7.556036e-05,
%      	m1 7.028195e-05, 
%     	oo1 7.824458e-05,&
%       rho1  6.531174e-05,
%      	ssa 0.003982e-04,&
%       m6 4.215566e-04,&
%       m8   5.620755e-04, 
%       mk3 2.134402e-04, 
%       s6 4.363323e-04,
%      	2sm2 1.503693e-04,&
%       2mk3  2.081166e-04,


