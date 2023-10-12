%--------------------------------------------------------------
%
%  Modify child grid topography such that it matches the interpolated
%  parent topography at the boundaries.
%
%  This script is for use with make_roms2roms.
%
%   (c) 2007 Jeroen Molemaker, UCLA
%--------------------------------------------------------------
%
if 0
  clear all
% ROMS parent and child grid directories
  pdir = '/avatar/nmolem/??/';
  cdir = '/avatar/nmolem/??/';
% ROMS parent and child grid files
  pgrid = '???_grd.nc';
  cgrid = '???_grd.nc';

% Only match to parent topography on open boundaries
  obcflag              = [1 1 1 1];      % open boundaries flag (1=open , [S E N W])

  pgrid = [pdir pgrid];
  cgrid = [cdir cgrid];
end
% End user-defined----------------------------------------------
%

%   Get topography data from childgrid
    lonc = ncread(cgrid,'lon_rho') ;%lonc = lonc-360;
    latc = ncread(cgrid,'lat_rho') ;
    hc   = ncread(cgrid,'hraw') ;
    maskc= ncread(cgrid,'mask_rho') ;
    [nxc,nyc]=size(hc);

  %  if parent_GLORYS

  % Get lat/lon from the child grid
      lonchd = double(ncread(cgrid,'lon_rho') );
      latchd = double(ncread(cgrid,'lat_rho') );

  %   Get lat/lon from bathy file
  %%% GLORYS longitude goes from -180 to 180
  %%% ROMS longitude goes from 0 to 360
      lonb = double(ncread(GLORYS_bathy,'longitude') );
      lonb(lonb < 0) = lonb(lonb < 0) + 360;
      latb = double(ncread(GLORYS_bathy,'latitude') );

      minlon = min(lonchd(:));
      maxlon = max(lonchd(:));
      minlat = min(latchd(:));
      maxlat = max(latchd(:));

      imin_tmp = find(lonb < minlon);
      imin_tmp2 = diff(imin_tmp);
      imin = find(imin_tmp2>1);

      imax_tmp = find(lonb > maxlon);
      imax_tmp2 = diff(imax_tmp);
      imax = find(imax_tmp2>1)+1;
      if isempty(imax)
        imax = imax_tmp(1);
      end

      jmin_tmp = find(latb < minlat);
      jmin = jmin_tmp(end);

      jmax_tmp = find(latb > maxlat);
      jmax = jmax_tmp(1);

      lon_tmp = lonb(imin:imax);
      lat_tmp = latb(jmin:jmax);

      [latp, lonp] = meshgrid(lat_tmp, lon_tmp);

      %lonp(lonp<0) = lonp(lonp < 0) + 360;
      hp   = double(ncread(GLORYS_bathy,'deptho') );
      hp = hp(imin:imax, jmin:jmax);
      maskp= double(ncread(GLORYS_bathy,'mask') );
      maskp= maskp(imin:imax, jmin:jmax, 1);
      [nxp,nyp] = size(hp);

 %   else
%   Get topography data from parent grid
 %     lonp = double(ncread(pgrid,'lon_rho') );
 %     latp = double(ncread(pgrid,'lat_rho') );
 %     hp   = double(ncread(pgrid,'h') );
 %     maskp= double(ncread(pgrid,'mask_rho') );
 %     lonp = lonp; % + 360;
 %     [nxp,nyp] = size(hp);
 %   end

%   Find minimal parent subgrid bounds
    lon0 = min(min(lonc));
    lon1 = max(max(lonc));
    lat0 = min(min(latc));
    lat1 = max(max(latc));

    spc = 2;
    g = lonp>=lon0&lonp<=lon1 & latp>=lat0&latp<=lat1;
    jmin = min(find(any(g ))) - spc;
    jmax = max(find(any(g ))) + spc;
    imin = min(find(any(g'))) - spc;
    imax = max(find(any(g'))) + spc;
    clear  g

    imin = max(imin,1  );
    jmin = max(jmin,1  );
    imax = min(imax,nxp);
    jmax = min(jmax,nyp);

   % for smode into westc grid
%   imin = imin + 199;
%   jmin = jmin + 180;
%   jmax = jmax - 187;


    hp    =   hp(imin:imax,jmin:jmax);
    lonp  = lonp(imin:imax,jmin:jmax);
    latp  = latp(imin:imax,jmin:jmax);
    maskp =maskp(imin:imax,jmin:jmax);
    [nxp,nyp]=size(hp);

    if 0
     plot(lonp,latp,'.k')
     hold on;plot(lonc,latc,'.r');hold off
     drawnow
     return
    end


%   Get interpolation coefficient to go to from (lonp,latp) to (lonc,latc).
    [elem,coef] = get_tri_coef(lonp,latp,lonc,latc,maskp);
    disp('get tri coef done')

%   Parent grid topo at child locations
    hpi = sum(coef.*hp(elem),3);
%   Parent mask at child locations
    maskpi = sum(coef.*maskp(elem),3);


  if 0
   dist = zeros(nxc,nyc,4);
   for i = 1:Mc   %% north south
    for j = 1:Lc    %% east west
     dist(i,j,1) =      i/Mc + (1-obcflag(1))*1e6; % South
     dist(i,j,2) = (Lc-j)/Lc + (1-obcflag(2))*1e6; % East
     dist(i,j,3) = (Mc-i)/Mc + (1-obcflag(3))*1e6; % North
     dist(i,j,4) =      j/Lc + (1-obcflag(4))*1e6; % West
    end
   end
   dist = min(dist,[],3);

   alpha = 0.5*tanh(100*(dist-0.03))+0.5; %% Feel free to play with this function.
   alpha = 0.5*tanh( 50*(dist-0.06))+0.5; %% Feel free to play with this function.
  else
   distance
  end

  %return

  hcn = alpha.*hc + (1-alpha).*hpi;
  maskcn = alpha.*maskc + (1-alpha).*maskpi;

  maskcn(maskcn>=0.5) = 1;
  maskcn(maskcn<0.5)  = 0;

  hcn(isnan(hcn)) = hmin;

%% Write modified topography to child grid
ncwrite(cgrid,'h',hcn);
ncwrite(cgrid,'mask_rho',maskcn);
ncwriteatt(cgrid,'h','Notes', ...
 'Topo has been modified to match the parent grid topo at the boundaries');
ncwriteatt(cgrid,'mask_rho','Notes', ...
 'Mask has been modified to match the parent grid Mask at the boundaries');


%% Visualize the modification
  sc0 = min(hcn(:));
  sc1 = max(max(hcn));
  subplot(2,2,1)
  pcolor(lonc,latc,hpi);caxis([sc0 sc1]);colorbar;shading flat
  title('Interpolated Parent Topo')
  subplot(2,2,2)
  pcolor(lonc,latc,hcn);caxis([sc0 sc1]);colorbar;shading flat
  title('Boundary Smoothed Child Topo')
  subplot(2,2,3)
  pcolor(lonc,latc,hcn-hpi);colorbar;shading flat
  title('Difference between Parent and child Topo');
  subplot(2,2,4)
  pcolor(lonc,latc,alpha);colorbar;shading flat
  title('Parent/Child transition function');
