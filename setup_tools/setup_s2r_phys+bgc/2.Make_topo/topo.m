
   gridfile = '/glade/scratch/bachman/ROMS_tools/Iceland0_BGC2/1.Make_grid/Iceland0_grd.nc';

%  Make
%  clear jet
%  colormap(jet(256))

   d2r = pi/180;
   r2d = 180/pi;
   r_earth = 6371315.;

   min_depth = 2;
%  1. Extract minimal set from topo file.  (plus 'width' size borders)
%  2. Establish bounding rectangle.
%  3. Based on weighting scale, determine maximum stencil.

%  TODO: add wrapper border whenever the grid stradles the dateline.
   % input parameters
for formask = 0:1

   if formask
     wd = 1.;  % for mask
%    climat = 'etopo';
%    climat = 'taidp200';
%    climat = 'taidp500';
     climat = 'srtm15';
%     climat = 'emodnet';
%    climat = 'srtm30';
%    climat = 'usgs_chesa';
   else
     wd=8;  % for hraw 8
%    climat = 'crm_west';
%    climat = 'srtm30';
     climat = 'srtm15';
%     climat = 'emodnet';
%    climat = 'etopo';
%    climat = 'palau';
%    climat = 'taidp200';
%    climat = 'taidp500';
%    climat = 'usgs_chesa';
   end
   wd = wd/2;


   global_grid = 1;
   if global_grid
    xg_glob = ncread(gridfile,'lon_rho');
    yg_glob = ncread(gridfile,'lat_rho');
    pn_glob = ncread(gridfile,'pn');
    h_glob = xg_glob + nan;
   end

   [nxg,nyg] = size(xg_glob);

   display(['Smoothing radius: ' num2str(wd) ' gridpoints'])

   if strcmp(lower(climat),'taidp500') & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from Taidp500m (smoothing radius ' num2str(wd) ')']);
   end
   if strcmp(lower(climat),'taidp200') & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from Taidp200m (smoothing radius ' num2str(wd) ')']);
   end
   if strcmp(lower(climat),'srtm30') & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from SRTM30 (Smoothing radius ' num2str(wd) ')']);
   end
   if strcmp(lower(climat),'srtm15') & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from SRTM15 V2.5 (Smoothing radius ' num2str(wd) ')']);
   end
   if strcmp(lower(climat),'emodnet') & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from EMODnet (Smoothing radius ' num2str(wd) ')']);
   end
   if strcmp(lower(climat),'etopo')  & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from ETOPO1 (Smoothing radius ' num2str(wd) ')']);
   end
   if strcmp(lower(climat),'usgs_chesa')  & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from USGS Coastal Relief 3-second data (Smoothing radius ' num2str(wd) ' )']);
   end
   if strcmp(lower(climat),'palau')  & ~formask
    ncwriteatt(gridfile,'hraw','source',['Raw bathymetry from Palau NOAA multi beam repository (Smoothing radius ' num2str(wd) ' )']);
   end


   nchx =  1;
   nchy =  1;
   chx = round([0:nchx]*nxg/nchx);
   chy = round([0:nchy]*nyg/nchy);
   ci = 1;
   cj = 1;

  xl0 = 1000;
  xl1 =-1000;
  yl0 = 1000;
  yl1 =-1000;
% for  ci = 3:3
%  for  cj = 2:2
  for  ci = 1:nchx
   for  cj = 1:nchy
     tic
     disp(['chunk ' num2str(ci) num2str(cj)])
     gj0 = chy(cj)+1;
     gj1 = chy(cj+1);
     gi0 = chx(ci)+1;
     gi1 = chx(ci+1);

%    ncg = netcdf(gridfile,'w');
     xg = xg_glob(gi0:gi1,gj0:gj1);
     yg = yg_glob(gi0:gi1,gj0:gj1);
     pn = pn_glob(gi0:gi1,gj0:gj1);
     [nxg,nyg] = size(xg);

     dg = 1./min(min(pn));
     width = wd*dg;  %% largest smoothing width for this patch

     %% figure out how large the patch of topo data needs to be
     xg_min = min(min(xg));
     xg_max = max(max(xg));
     yg_min = min(min(yg));
     yg_max = max(max(yg));

     coslat = cos(d2r*max(max(abs(yg))));
     dellon = r2d*width/(coslat*r_earth);
     dellat = r2d*width/(r_earth);

     xg_min = xg_min-1.2*dellon*2;
     xg_max = xg_max+1.2*dellon*2;
     yg_min = yg_min-1.2*dellat*2;
     yg_max = yg_max+1.2*dellat*2;

     xl0 = min(xl0,xg_min);   %% Just for plotting
     xl1 = max(xl1,xg_max);
     yl0 = min(yl0,yg_min);
     yl1 = max(yl1,yg_max);

     switch lower(climat);
      case{'taidp200'}
        disp('Using Taidp200')
        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/Taidp200.nc'; %% from
        topo_lon = ncread(datafile,'longitude');
        topo_lat = ncread(datafile,'latitude');
      case{'taidp500'}
        disp('Using Taidp500')
        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/Taidp500.nc'; %% from
        topo_lon = ncread(datafile,'longitude');
        topo_lat = ncread(datafile,'latitude');
      case{'crm_west'}
        disp('Using CRM_west')
        datafile = '/batavia/nmolem/OBSERV/TOPO/crm_west.nc'; %% from
        topo_lon = ncread(datafile,'x');
        topo_lat = ncread(datafile,'y');
       %topo_lon = topo_lon + 360;
      case{'usgs_chesa'}
        disp('Using USGS data set')
        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/usgs_chesapeake.nc'; %% from -77 to -74 degrees
        topo_lon = ncread(datafile,'lon');
        topo_lat = ncread(datafile,'lat');
      case{'etopo'}
        disp('Using Etopo1 data set')
%       nct = netcdf('/batavia/nmolem/OBSERV/TOPO/ETOPO2v2c_f4.nc'); %% from -180 180 degrees
%       nct = netcdf('/batavia/nmolem/OBSERV/TOPO/etopo2_0_360.nc'); %% from 0 360 degrees
%       topo_lon = nct{'x'}(:);
%       topo_lat = nct{'y'}(:);

        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/ETOPO1_180.nc'; %% from -180 180 degrees
        topo_lon = ncread(datafile,'lon');
        topo_lat = ncread(datafile,'lat');
      case{'gebco'}
        disp('Using Gebco data set')
        nct = netcdf('/batavia/nmolem/OBSERV/TOPO/topo_gebco_GRANCAN.nc');
        topo_lon = nct{'lon'}(:);
        topo_lat = nct{'lat'}(:);
      case{'srtm30'}
        disp('Using SCRIPS data set')
%       datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/SRTM30_180.nc';
        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/SRTM30_0_360_new.nc';
        topo_lon = ncread(datafile,'longitude');%   - 360;
        topo_lat = ncread(datafile,'latitude');
      case{'palau'}
        datafile = '/avatar/nmolem/batavia/nmolem/OBSERV/TOPO/palau_bathy.nc';
        topo_lon = ncread(datafile,'lon');
        topo_lat = ncread(datafile,'lat');
	topo_lon = topo_lon(:,1);
	topo_lat = topo_lat(1,:);
      case{'srtm15'}
        disp('Using the SRTM15 data set')
        datafile = '/glade/scratch/bachman/ROMS_tools/DATASETS/SRTM15_V2.5.nc';
        topo_lon = ncread(datafile,'lon');   %%% adding 360 for locations west of dateline
        topo_lon(topo_lon < 0) = topo_lon(topo_lon < 0) + 360;
        topo_lat = ncread(datafile,'lat');
        tmp = topo_lon(1:43200);
        topo_lon(1:43200) = topo_lon(43201:86400);
        topo_lon(43201:86400) = tmp;
      case{'emodnet'}
        disp('Using the EMODnet data set')
        datafile = '/glade/scratch/bachman/ROMS_tools/code/DATASETS/EMODnet_C2.nc';
        topo_lon = ncread(datafile,'lon');   %%% adding 360 for locations west of dateline
        topo_lon(topo_lon < 0) = topo_lon(topo_lon < 0) + 360;
        topo_lat = ncread(datafile,'lat');
     end

     dx_topo = r_earth*(topo_lon(2)-topo_lon(1))*d2r;
     dy_topo = r_earth*(topo_lat(2)-topo_lat(1))*d2r;

%     [min(topo_lon) xg_min]
%     [max(topo_lon) xg_max]

     max(topo_lat)
     min(topo_lat)
     max(topo_lon)
     min(topo_lon)
     yg_min
     yg_max
     xg_min
     xg_max

     if max(topo_lat)<yg_min|min(topo_lat)>yg_max
       disp('chunk not in range')
       continue
     end
     if max(topo_lon)<xg_min|min(topo_lon)>xg_max
       disp('chunk not in range')
       continue
     end
     if ~strcmp(lower(climat),'sands')
      ibt = find( topo_lon < xg_min,100000,'first')-1;
      ibt = ibt(end);
      iet = find( topo_lon > xg_max,1,'first')+1;
      jbt = find( topo_lat > yg_min,1,'first')-1;
      jet = find( topo_lat < yg_max,1,'last')+1;
      ibt = max(ibt,1); iet = min(iet,length(topo_lon));
      jbt = max(jbt,1); jet = min(jet,length(topo_lat));
%     if yg_min<min(topo_lat)+0.04
%      jbt = 1;
%      yg(yg<min(topo_lat)+0.04) = min(topo_lat)+ 0.04;
%     end
%     if xg_min<min(topo_lat)+0.04
%      ibt = 1;
%      xg(xg<min(topo_lon)+0.04) = min(topo_lon)+ 0.04;
%     end

      topo_lon = topo_lon(ibt:iet);
      topo_lat = topo_lat(jbt:jet);
     end



  switch lower(climat);
    case{'taidp200'}
      hin = -ncread(datafile,'elevation',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'taidp500'}
      hin = -ncread(datafile,'elevation',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'crm_west'}
%     hin = nct{'topo'}(jbt:jet,ibt:iet);
      hin = -ncread(datafile,'z',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'usgs_chesa'}
%     hin = nct{'topo'}(jbt:jet,ibt:iet);
      hin = ncread(datafile,'topo',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'etopo'}
      hin = -ncread(datafile,'z',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'gebco'}
      hin = -nct{'topo'}(jbt:jet,ibt:iet);
    case{'srtm30'}
%     hin = -nct{'elevation'}(jbt:jet,ibt:iet);
      hin = -ncread(datafile,'elevation',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
%     hin(hin>7500) = 7500;disp('limiting srtm30 topo at 7500m');
    case{'sands'}
      [hin,topo_lat,topo_lon] = extract_1m([yg_min-0.02 yg_max+0.02 xg_min-0.02 xg_max+0.02],1);
       topo_lat = topo_lat(1:end-1);
       hin = -hin(2:end,:);
    case{'palau'}
      hin = ncread(datafile,'fdept',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
    case{'srtm15'}
      if ibt<43200
        hin = -ncread(datafile,'z',[(ibt+43200) jbt],[iet-ibt+1 jet-jbt+1]);
      else
        hin = -ncread(datafile,'z',[(ibt-43199) jbt],[iet-ibt+1 jet-jbt+1]);
      end
    case{'emodnet'}
      hin = -ncread(datafile,'elevation',[ibt jbt],[iet-ibt+1 jet-jbt+1]);
   end
   [nxt,nyt] = size(hin);

   if ~formask
%   hin(hin<5) = 5;
%   hin(hin<5) = 5;
   end

   if sum(~isnan(hin(:)))<10
     display('no data for this chunk')
     continue
   end

   if 0  %% testing
     [xt2,yt2] = meshgrid(xt,yt);
     hin = xt2;
   end

   ill = xg;
   jll = yg;
   lon_min = min(topo_lon);
   lat_min = min(topo_lat);
   for i = 1: nxg  %% Locate a index in the topo file that is close to the grid point
    for j = 1: nyg
      if lat_min>yg(i,j)
        jll(i,j) = 1;
      else
        jll(i,j) = find( topo_lat < yg(i,j),1,'last');
      end
      if lon_min>xg(i,j)
        ill(i,j) = 1;
      else
        ill(i,j) = find( topo_lon < xg(i,j),1,'last');
      end
    end
   end
   sl = size(topo_lat)


   %% minimal distance of dx is:
   %%  xt and yt should be monotonous and increasing
    lcoslat = cos(d2r*max(abs(topo_lat)));

    isz = 2*ceil(width/dx_topo/lcoslat)
    jsz = 2*ceil(width/dy_topo)
%   siz = max(isz,jsz);
%   isz = siz;jsz=siz;

 if isz<3 & jsz<3
   display('Doing linear interpolation');
   [xt,yt] = meshgrid(topo_lon,topo_lat);
%  xt = xt';
%  yt = yt';
%  hraw = interp2(xt,yt,hin,xg,yg,'spline');
   hraw = interp2(xt,yt,hin',xg,yg,'linear');
 else
   display('Doing Weighted average');

   shfti = isz/2 -1;
   shftj = jsz/2 -1;
   weight = zeros(isz,jsz);
   hraw   = zeros(nxg,nyg) + nan;
   for i = 1:nxg
    for j = 1:nyg
      ib = ill(i,j)      -shfti;
      ie = ill(i,j)+isz-1-shfti;

      jb = jll(i,j)      -shftj;
      je = jll(i,j)+jsz-1-shftj;

      ib = max(ib,1);ie = min(ie,nxt);
      jb = max(jb,1);je = min(je,nyt);

      xtl = topo_lon(ib:ie);
      ytl = topo_lat(jb:je);
      [xtl,ytl] = meshgrid(xtl,ytl);
      xtl = xtl'; ytl = ytl';
      r = gc_dist(xg(i,j),yg(i,j),xtl,ytl,1) ; % the 1 argument is a radian to degree conversion
      lwidth = wd./pn(i,j);  %% local cell size
%     [width lwidth]
      w = (1-(r/lwidth).^2).^2;
      w(r>lwidth) = 0;

      hlocal = hin(ib:ie,jb:je);

%     if 0 % used for crm
%     w(isnan(hlocal)) = 0;
%     hlocal(isnan(hlocal)) = 0;
%     ratio = sum(w(:)==0)/length(w(:));
%     if ratio>0.5
%       continue
%     end
%     end

      w = w/sum(sum(w));
%     return
%     if isnan(w)>0
%     end

      if min([min(w(1,:)) min(w(end,:)) min(w(:,1)) min(w(:,end))]) >0
        figure(3);contour(w,[0 0],'k');drawnow
        error 'topo sampling region too small'
      end
      hraw(i,j) = sum(w(:).*hlocal(:));


    end
    if mod(i,10)==0
      hraw(hraw<-0.2) = -0.2;
      figure(1); mypcolor(xg,yg,hraw);
      pause(0.01)
    end

   end
 end

% hraw(hraw<-0.2) = -0.2;
   %hraw(hraw<min_depth) = min_depth;
% hraw(hraw==0) = nan;
   %hraw(isnan(hraw)) = 0;

% hraw(hraw<0) = 0;
  mask = 1 + 0.*hraw;
  mask(hraw<=0.11) = 0;

  figure(2); mypcolor(xg,yg,mask);
  pause(0.01)

%   if formask
     disp('writing mask')
     ncwrite(gridfile,'mask_rho',mask,[gi0 gj0]);
%    ncg{'mask_rho'}(gi0:gi1,gj0:gj1)= mask;
%   else
     disp('writing topo')
     ncwrite(gridfile,'hraw',hraw,[gi0 gj0]);
%   end

%  if global_grid
%   h_glob(gi0:gi1,gj0:gj1) =  hraw;
%   figure(2);
%   mypcolor(xg_glob,yg_glob,h_glob);xlim([xl0 xl1]);ylim([yl0 yl1])
%  end
 toc
end  %% chj loop
end  %% chi loop


end % if formask
