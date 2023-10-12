% correction_on_frcbgc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  2019, Pierre Damien (UCLA)
%
%
% Create BGC file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear all
 close all

grdname  = '/glade/scratch/bachman/ROMS_tools/setup_r2r_phys+bgc/1.Make_grid/Iceland1_grd.nc';
root_name= '/glade/scratch/bachman/ROMS_tools/setup_r2r_phys+bgc/1.Make_grid/Iceland1';

coarse_frc = 0;

   frcname = [root_name '_frc_bgc.nc']
   if exist(frcname)
     delete(frcname)
   end

   disp(' Read in the target grid')

   if coarse_frc
     lon = ncread(grdname,'lon_coarse');
     lat = ncread(grdname,'lat_coarse');
   else
     lon = ncread(grdname,'lon_rho');
     lat = ncread(grdname,'lat_rho');
   end
   [nx,ny] = size(lon);

   % extend data periodically in longitude
   % if there are grid locations out of range

   grd.lon = lon;
   grd.lat = lat;

   lon0 = min(lon(:));
   lon1 = max(lon(:));
   lat0 = min(lat(:));
   lat1 = max(lat(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp(' ')
   disp(' Read in the pco2 file ')

   datname = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/spco2_aco2_clim_ETH_SOM-FFN_CDIAC_G05_corr.nc';
   lon_frc = ncread(datname,'lon'); lon_frc = lon_frc(:,1) ;
   lat_frc = ncread(datname,'lat'); lat_frc = lat_frc(1,:) ;
   dat.nx = length(lon_frc);
   dat.ny = length(lat_frc);

   % figure out periodic extention of data
   if lon0<min(lon_frc)
     disp('extending west')
     ext_west = 1;
     i0 = find(lon_frc-360<lon0,1,'last');
   else
     ext_west = 0;
     i0 = find(lon_frc<lon0,1,'last');
   end
   if lon1>max(lon_frc)
     disp('extending east')
     ext_east = 1;
     i1 = find(lon_frc+360>lon1,1,'first');
   else
     ext_east = 0;
     i1 = find(lon_frc>lon1,1,'first');
   end

   j0 = find(lat_frc<lat0,1,'last');
   j1 = find(lat_frc>lat1,1,'first');
   fny = j1-j0+1;

   if ext_west
     lon_frc= [lon_frc(i0:end)'-360 lon_frc(1:i1)'];
   elseif ext_east
     lon_frc= [lon_frc(i0:end)' lon_frc(1:i1)'+360];
   else
     lon_frc= lon_frc(i0:i1);
   end
   lat_frc = lat_frc(j0:j1);

   create_frc_BGC(grdname,frcname,coarse_frc,'pco2');
   pco2_time = [0.5:11.5]*365.25/12; % Months of equal length
   ncwrite(frcname,'pco2_time',pco2_time);
      % Loop over data files
   for i = 1:12
     disp([' Processing month : ' num2str(i)])
       if ext_west | ext_east
       fnx1 = dat.nx-i0+1;
       fnx2 = i1;
       data1 = ncread(datname,'pco2_air',[i0 j0 i],[fnx1 fny 1]);
       data2 = ncread(datname,'pco2_air',[1  j0 i],[fnx2 fny 1]);
       data = [data1' data2']';
     else
       fnx = i1-i0+1;
       data = ncread(datname,'pco2_air',[i0 j0 i],[fnx fny 1]);
     end
     data = inpaint_nans(data,2);
     co2 = interp2(lon_frc,lat_frc,data',lon,lat,'makima');
     ncwrite(frcname,'pco2_air',co2,[1 1 i])
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp(' ')
   disp(' Read in the iron file ')

   datname = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/solFe_scenario4_current_gx1v6.nc';
   factor = 1.79*10^6; % conversion from dust [kg/(m2 s)] to iron [nmol/(cm2 s)]:
   lon_frc = ncread(datname,'TLONG'); %lon_frc(lon_frc>0) = lon_frc(lon_frc>0)-360 ;
   lat_frc = ncread(datname,'TLAT');

   create_frc_BGC(grdname,frcname,coarse_frc,'iron');
   iron_time = [0.5:11.5]*365.25/12; % Months of equal length
   ncwrite(frcname,'iron_time',iron_time);
      % Loop over data files
   for i = 1:12
     disp([' Processing month : ' num2str(i)])
     data = ncread(datname,'DSTSF',[1 1 i],[inf inf 1]);
     data(data==0) = NaN ;
     iron = griddata(double(lon_frc),double(lat_frc),double(data),lon,lat);
     iron = inpaint_nans(iron,2);
     iron = iron * factor ;
     ncwrite(frcname,'iron',iron,[1 1 i])
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp(' ')
   disp(' Read in the dust file ')

   datname = '/glade/scratch/bachman/ROMS_tools/DATASETS/BGC/dst79gnx_gx1v6_090416.nc';
   factor = 1 ;
   lon_frc = ncread(datname,'TLONG');% lon_frc(lon_frc>0) = lon_frc(lon_frc>0)-360 ;
   lat_frc = ncread(datname,'TLAT');

   create_frc_BGC(grdname,frcname,coarse_frc,'dust');
   dust_time = [0.5:11.5]*365.25/12; % Months of equal length
   ncwrite(frcname,'dust_time',dust_time);
      % Loop over data files
   for i = 1:12
     disp([' Processing month : ' num2str(i)])
     data = ncread(datname,'DSTSF',[1 1 i],[inf inf 1]);
     data(data<0) = NaN ;
     dust = griddata(double(lon_frc),double(lat_frc),double(data),lon,lat);
     dust = inpaint_nans(dust,2);
     dust = dust * factor ;
     ncwrite(frcname,'dust',dust,[1 1 i])
   end


correction_on_FRCtime;






