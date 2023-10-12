% fill_frc_frc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make climatologic SSS field forcing file using WOA data
%
%  2022, Jeroen Molemaker, Pierre Damien, UCLA
%
%
%%%%%%%%%%%%%%%%%%%%% USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
%
%  frc climatology file names:

frc_dir = '/paracas/DATASETS/WOA18/salinity/';

grdname  = '/paracas/nmolem/PACHUG/pachug_grd.nc';
root_name= '/paracas/nmolem/PACHUG/pachug';

grdname  = '/paracas/nmolem/NWPAC/nwpac_grd.nc';
root_name= '/paracas/nmolem/NWPAC/nwpac';

grdname  = '/paracas/nmolem/NEPAC/nepac_grd.nc';
root_name= '/paracas/nmolem/NEPAC/nepac';

grdname  = '/paracas/nmolem/LUZON/luzon_grd.nc';
root_name= '/paracas/nmolem/LUZON/luzon';

grdname  = '/paracas/nmolem/NORMAR/normar_grd.nc';
root_name= '/paracas/nmolem/NORMAR/normar';

grdname  = '/paracas/nmolem/GREEN/green_grd.nc';
root_name= '/paracas/nmolem/GREEN/green';

coarse_frc = 1;

%
%%%%%%%%%%%%%%%%%%% END USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%

%
   ssslist = dir([frc_dir 'woa18_decav_s*']);
   nfiles = length(ssslist);

   % find the right files

   % Read in data coordinates and set data trim.
   disp(' ')
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

   disp(' ')
   disp(' Read in the data grid')

   datname = [frc_dir ssslist(1).name];
   lon_frc = ncread(datname,'lon');
   lat_frc = ncread(datname,'lat');
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


   if nfiles~=12
    disp('not 12 monthly files')
    return
   end

   frcname = [root_name '_frc_sss.nc'] 
   if exist(frcname)
     delete(frcname)
   end
   create_frc_sss(grdname,frcname,coarse_frc);

   sss_time = [0.5:11.5]*365.25/12; % Mornths of equal length
   ncwrite(frcname,'sss_time',sss_time);

   % Loop over data files
   for i = 1:nfiles
     datname = [frc_dir ssslist(i).name];
     disp([' Processing: ' ssslist(i).name])

     if ext_west | ext_east
       fnx1 = dat.nx-i0+1;
       fnx2 = i1;
       sss1 = ncread(datname,'s_an',[i0 j0 1 1],[fnx1 fny 1 1]); 
       sss2 = ncread(datname,'s_an',[1  j0 1 1],[fnx2 fny 1 1]); 
       data.sss = [sss1' sss2']';
     else
       fnx = i1-i0+1;
       data.sss = ncread(datname,'s_an',[i0 j0 1 1],[fnx fny 1 1]); 
     end

     data.sss = inpaint_nans(data.sss,2);
     sss = interp2(lon_frc,lat_frc,data.sss',lon,lat,'makima');
     ncwrite(frcname,'sss',sss,[1 1 i])

   end
   return


