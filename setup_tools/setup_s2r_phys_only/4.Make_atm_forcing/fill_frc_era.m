% fill_frc_frc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make a ROMS bulk forcing file using ERA5 hourly data
%
%  2022, Jeroen Molemaker, Pierre Damien, UCLA
%
%  Possible future work
% - Extend lon_frc (and data) when grids straddle the prime meridian
% - Investigate over which time period the fluxes are integrated (centered or shifted)
% - Add the effect of surface pressure on Qair and other things
% - Maybe force the model with Tair and Dew point instead of Humidity
%
%%%%%%%%%%%%%%%%%%%%% USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
%
%  frc climatology file names:

frc_dir = '/glade/scratch/bachman/ERA5/NA/';
%frc_dir = '/glade/scratch/bachman/ERA5/CT/';

% Set a date range for the forcing file
start_date = datenum(2012,1,2);
end_date   = datenum(2012,1,12);
%start_date = datenum(2019,5,13);
%end_date   = datenum(2019,5,15);

%grdname = '/paracas/nmolem/PACHUG/pachug_grd.nc';
%disname = '/paracas/nmolem/PACHUG/pachug_cdist.mat';
%rivname = '/paracas/nmolem/PACHUG/pachug_swf.nc';
%root_name='/paracas/nmolem/PACHUG/pachug';

%grdname  = '/paracas/nmolem/NWPAC/nwpac_grd.nc';
%disname  = '/paracas/nmolem/NWPAC/nwpac_cdist.mat';
%root_name= '/paracas/nmolem/NWPAC/nwpac';

%grdname   ='/paracas/nmolem/NORMAR/normar_grd.nc';
%disname   ='/paracas/nmolem/NORMAR/normar_cdist.mat';
%root_name ='/paracas/nmolem/NORMAR/normar';

%grdname   ='/paracas/nmolem/GREEN/green_grd.nc';
%disname   ='/paracas/nmolem/GREEN/green_cdist.mat';
%root_name ='/paracas/nmolem/GREEN/green';

%grdname  = '/paracas/nmolem/NEPAC/nepac_grd.nc';
%disname  = '/paracas/nmolem/NEPAC/nepac_cdist.mat';
%root_name='/paracas/nmolem/NEPAC/nepac';

%grdname  = '/paracas/nmolem/LUZON/luzon_grd.nc';
%disname  = '/paracas/nmolem/LUZON/luzon_cdist.mat';
%root_name= '/paracas/nmolem/LUZON/luzon';

%grdname  = '/paracas/nmolem/ANGAUR/angaur_grd.nc';
%disname  = '/paracas/nmolem/ANGAUR/angaur_cdist.mat';
%root_name= '/paracas/nmolem/ANGAUR/angaur';

grdname= '/glade/scratch/bachman/ROMS_tools/setup_s2r_phys_only/1.Make_grid/Wales0_grd.nc';
root_name= '/glade/scratch/bachman/ROMS_tools/setup_s2r_phys_only/1.Make_grid/Wales0';

disname  = '/paracas/nmolem/WESTC/westc_cdist.mat';

add_coarse_grd(grdname);

swcorrname = [frc_dir 'SSR_correction.nc'];


coarse_frc   = 1; % forcing files at half the resolution of the grid
wind_dropoff = 0;
add_rivers   = 0; % Adds river runoff as additional precipitation

%
%%%%%%%%%%%%%%%%%%% END USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%

%
   eralist = dir([frc_dir 'ERA5*']);
   nfiles = length(eralist);

   start_date
   end_date
   % find the right files
   stimes = zeros(nfiles,1);
   etimes = zeros(nfiles,1);   %%%%
   for i=1:nfiles
     datname = [frc_dir eralist(i).name];
     times = double(ncread(datname,'time'))/24 + datenum(1900,1,1)
     stime(i) = times(1)
     etime(i) = times(end)
   end

   t0 = find(stime<start_date,1,'last');
   t1 = find(etime>end_date,1,'first');   %% stime

   %t0 = t0+1
   % trim list of filenames
   eralist = eralist(t0:t1);
   nfiles = length(eralist);
%  eralist.name
%  return

   % Read in data coordinates and set data trim.
   disp(' ')
   disp(' Read in the target grid')
   if coarse_frc
%    if ~exist(grdname,'file')
       disp('adding coarse lon/lat to grid')
       add_coarse_grd(grdname);
%    end
     disp('Reading coarse grid')
     lon = ncread(grdname,'lon_coarse');
     lat = ncread(grdname,'lat_coarse');
     ang = ncread(grdname,'angle_coarse');
   else
     disp('Reading full grid')
     lon = ncread(grdname,'lon_rho');
     lat = ncread(grdname,'lat_rho');
     ang = ncread(grdname,'angle');
   end
   [nx,ny] = size(lon);

   lon = mod(lon,360); %% ERA5 forcing data is in the 0 to 360 longitude range
   cosa = cos(ang);
   sina = sin(ang);

   grd.lon = lon;
   grd.lat = lat;

   lon0 = min(lon(:))
   lon1 = max(lon(:))
   lat0 = min(lat(:))
   lat1 = max(lat(:))

   disp(' ')
   disp(' Read in the data grid')

   datname = [frc_dir eralist(1).name];

   lon_frc = ncread(datname,'longitude');
   lon_frc(lon_frc < 0) = lon_frc(lon_frc < 0) + 360;   % ERA5 forcing data is in the 0 to 360 longitude range
   lat_frc = ncread(datname,'latitude');

   i0 = find(lon_frc<lon0,1,'last');
   i1 = find(lon_frc>lon1,1,'first');

   % ERA5 is written in upside down order, latitude is decreasing

   j0 = find(lat_frc>lat1,1,'last');
   j1 = find(lat_frc<lat0,1,'first');
   fnx = i1-i0+1;
   fny = j1-j0+1;

   lon_frc = lon_frc(i0:i1);
   lat_frc = flipud(lat_frc(j0:j1));

   sst = ncread(datname,'sst',[i0 j0 1],[fnx fny 1]);
   mask = 1+0*sst;
   mask(isnan(sst)) = 0;
   mask = fliplr(mask); % to deal with upside down ERA5 data


   if wind_dropoff
     if ~exist(disname)
       cdist = comp_cdist(grdname,disname,coarse_frc);
     else
       load(disname)
     end
     cdist = cdist/1e3;
%    mult = 1-0.6*exp(-0.08*cdist);
     mult = 1-0.4*exp(-0.08*cdist);
   end

   % prepare for short wave radiation correction
   corr_time = ncread(swcorrname,'time');
   swr_mult = zeros(nx,ny,12);
   for i=1:12
     corr_swr  = ncread(swcorrname,'ssr_corr',[i0 j0 i],[fnx fny 1]);
     corr_swr = fliplr(corr_swr);
     corr_swr(mask<1) = nan;
     corr_swr = inpaint_nans(corr_swr,2);
     swr_mult(:,:,i) = interp2(lon_frc,lat_frc,corr_swr',lon,lat,'makima');
   end

   if add_rivers
     % prepare for the addition of river runoff to rain data
     riv_time = ncread(rivname,'swf_time');
     riv_flux = ncread(rivname,'swf_flux');
   end



   % Loop over data files
   for i = 1:nfiles
     datname = [frc_dir eralist(i).name];
     disp([' Processing: ' eralist(i).name])

     dat_time = ncread(datname,'time');
     nrecord = length(dat_time);

     date_num = double(dat_time(1)/24) + datenum(1900,1,1);
     label = datestr(date_num,'YYYYmm')

     frcname = [root_name '_frc.' label '.nc']

     data.datname = datname;
     data.lon = lon_frc;
     data.lat = lat_frc;
     data.i0 = i0;
     data.j0 = j0;
     data.fnx = fnx;
     data.fny = fny;
     data.mask = mask;

     if exist(frcname)
       delete(frcname)
     end
     create_frc_bulk(grdname,frcname,coarse_frc);

     for irec = 1:nrecord
        disp(['Record: ' num2str(irec)])

	% ---- time -----
	time = ncread(datname,'time',[irec],[1]);

	% translate to days since 2000,1,1
%       days = double(dat_time(1)/24) + datenum(1900,1,1) - datenum(2000,1,1)
        days = double(time)/24. + datenum(1900,1,1) - datenum(2000,1,1);
	ncwrite(frcname,'time',days,[irec]);

	% ---- 10 meter winds -----
	u = get_frc_era(data,grd,'u10',irec,'makima');
	v = get_frc_era(data,grd,'v10',irec,'makima');

	% rotate to grid angles
        ugrid = cosa.*u + sina.*v;
        vgrid =-sina.*u + cosa.*v;

	if wind_dropoff
          ugrid = ugrid.*mult;
          vgrid = vgrid.*mult;
        end

	ncwrite(frcname,'uwnd',ugrid,[1 1 irec]);
	ncwrite(frcname,'vwnd',vgrid,[1 1 irec]);

	% ---- Incoming Radiation -----
	swr = get_frc_era(data,grd,'ssr',irec,'linear');  % downward_shortwave_flux [J/m2]
	lwr = get_frc_era(data,grd,'strd',irec,'linear'); % downward_longwave_flux [J/m2]

	% Translate to fluxes. ERA5 stores values integrated over 1 hour
	% Are these centered around the current time?
	swr = swr/3600;
	lwr = lwr/3600;

	% ---- Correction to swr using  COREv2 dataset -----
	% This is a multiplicative correcting, trying to account for
	% for errors in the ERA5 cloud cover
	% temporal interpolation in a monthly climatology
	yr_day = mod(days,365.25);
	if yr_day>=corr_time(12)
           id0 = 12;
           id1 =  1;
	   cf  = (yr_day - corr_time(12))/(corr_time(1)+365.25-corr_time(12));
        elseif yr_day<=corr_time(1)
           id0 = 12;
           id1 =  1;
	   cf  = (yr_day - corr_time(12)+365.25)/(corr_time(1)+365.25-corr_time(12));
        else
           id0 = find(yr_day>corr_time,1,'last');
	   id1 = id0+1;
	   cf  = (yr_day - corr_time(id0))/(corr_time(id1)-corr_time(id0));
	end

	swrd= swr;
	swr_cr = cf*swr_mult(:,:,id0) + (1-cf)*swr_mult(:,:,id1);
	swr= swr.*swr_cr;

	ncwrite(frcname,'swrad',swr,[1 1 irec]);
	ncwrite(frcname,'lwrad',lwr,[1 1 irec]);

	% ---- Absolute humidity -----

	t2m = get_frc_era(data,grd,'t2m',irec,'linear'); % 2 meter air temp [K]
	d2m = get_frc_era(data,grd,'d2m',irec,'linear'); % 2 meter air dew point temp [K]
	t2m = t2m - 273.15; % K to C
	d2m = d2m - 273.15; % K to C
        Qair=(exp((17.625*d2m)./(243.04+d2m))./exp((17.625*t2m)./(243.04+t2m)));  % Relative humidity fraction

        % Relative to absolute humidity assuming constand pressure

        patm=1010.0;

        cff=(1.0007+3.46e-6*patm).*6.1121.*exp(17.502*t2m./(240.97+t2m));
        cff=cff.* Qair;
        qair =0.62197.*(cff./(patm-0.378.*cff)) ;

	ncwrite(frcname,'qair',qair,[1 1 irec]);
	ncwrite(frcname,'Tair',t2m, [1 1 irec]);

	% ---- Rain ------
	rain = get_frc_era(data,grd,'tp',irec,'linear'); % precipitation [m]

	% Translate to fluxes. ERA5 stores values integrated over 1 hour
	% Are these centered around the current time?
	if 0 % output in m/s
          rain = rain/3600;
        else % output in cm/day
	  rain = rain*100*24;
        end

	% ---- Correction to rain using Dai and Trenberth data  -----
	% Temporal interpolation in a monthly climatology
        if add_rivers
          yr_day = mod(days,365.25);
          if yr_day>=riv_time(12)
            id0 = 12;
            id1 =  1;
	    cf  = (yr_day -  riv_time(12))/( riv_time(1)+365.25- riv_time(12));
          elseif yr_day<= riv_time(1)
            id0 = 12;
            id1 =  1;
	   cf  = (yr_day -  riv_time(12)+365.25)/( riv_time(1)+365.25- riv_time(12));
          else
            id0 = find(yr_day>riv_time,1,'last');
	    id1 = id0+1;
	    cf  = (yr_day -  riv_time(id0))/( riv_time(id1)- riv_time(id0));
	  end
          rain = rain + cf*riv_flux(:,:,id0) + (1-cf)*riv_flux(:,:,id1);
        end

	ncwrite(frcname,'rain',rain, [1 1 irec]);

     end
   end %


