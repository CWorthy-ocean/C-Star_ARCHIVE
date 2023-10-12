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
%

%   Get topography data from childgrid
    lonc = ncread(cgrid,'lon_rho') ;%lonc = lonc-360;
    latc = ncread(cgrid,'lat_rho') ;
    hc   = ncread(cgrid,'h') ;
    maskc= ncread(cgrid,'mask_rho') ;
    [nxc,nyc]=size(hc);

%   Get topography data from parent grid
    lonp = double(ncread(pgrid,'lon_rho') );
    latp = double(ncread(pgrid,'lat_rho') );
    hp   = double(ncread(pgrid,'h') );
    maskp= double(ncread(pgrid,'mask_rho') );
    lonp = lonp; % + 360;
    [nxp,nyp] = size(hp);

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
