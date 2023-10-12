 function [lon4,lat4,pm,pn,ang,lone,late] = easy_grid(nx,ny,size_x,size_y,tra_lon,tra_lat,rot);
 %
 %   Easy grid makes an rectangular, orthogonal grid with minimal gridsize variation
 %   It uses an Mercator projection around the equator and then rotates the sphere
 %   around its 3 axis to position the grid where it is desired.
 %   
 %   Inputs:
 %   nx:      Number of grid point in the x direction
 %   ny:      Number of grid point in the y direction
 %   size_x:  Domain size in x-direction
 %   size_y:  Domain size in y-direction
 %   tra_lon: Desired longitude of grid center
 %   tra_lat: Desired latitude of grid center
 %   rot:     Rotation of grid direction (0: x direction is west-east)
 %
 %   Example:  > [lon,lat] = easy_grid(30,20,4e6,3e6,-180,59,0)
 %
 %  (c) 2008, Jeroen Molemaker, UCLA


   r_earth = 6371315.;

   %% Mercator projection around the equator

   if (size_y>size_x) 
     length = size_y; nl = ny;
     width  = size_x; nw = nx;
   else
     length = size_x; nl = nx;
     width  = size_y; nw = ny;
   end
   
   dlon = length/r_earth;
%  lon1d = dlon*[0:1:nl]/nl - dlon/2;  
   lon1d = dlon*[-0.5:1:nl+0.5]/nl - dlon/2;  
   mul = 1;
   dlat = width/r_earth;
   for it = 1:100
    y1 = log(tan(pi/4-dlat/4));
    y2 = log(tan(pi/4+dlat/4));
    y = (y2-y1)*[-0.5:1:nw+0.5]/nw + y1;  
    lat1d = 2*atan(exp(y)) - pi/2; 
    lat1d = atan(sinh(y));
    dlat_cen = 0.5*(lat1d(round(nw/2)+1)-lat1d(round(nw/2)-1));
    dlon_cen = dlon/nl;
    mul = dlat_cen/dlon_cen*length/width*nw/nl;
%   lat1d = lat1d/mul;
    dlat = dlat/mul;
   end

   
   lon1de= dlon*[-1:1:nl+1]/nl - dlon/2;  
   ye= (y2-y1)*[-1:1:nw+1]/nw + y1;  
   lat1de= 2*atan(exp(ye)) - pi/2; 
   lat1de = atan(sinh(ye));
   lat1de= lat1de/mul;

   [lon1,lat1] = meshgrid(lon1d,lat1d);
   [lone,late] = meshgrid(lon1de,lat1de);
   lonu = 0.5*(lon1(:,1:end-1)+lon1(:,2:end));
   latu = 0.5*(lat1(:,1:end-1)+lat1(:,2:end));
   lonv = 0.5*(lon1(1:end-1,:)+lon1(2:end,:));
   latv = 0.5*(lat1(1:end-1,:)+lat1(2:end,:));

   if (size_y>size_x) 
    [lon1,lat1] = rot_sphere(lon1,lat1,90);
    [lonu,latu] = rot_sphere(lonu,latu,90);
    [lonv,latv] = rot_sphere(lonv,latv,90);
    [lone,late] = rot_sphere(lone,late,90);
    lon1 = flipdim(lon1,1)';
    lat1 = flipdim(lat1,1)';
    lone = flipdim(lone,1)';
    late = flipdim(late,1)';

    lonu_tmp= flipdim(lonv,1)';
    latu_tmp = flipdim(latv,1)';
    lonv = flipdim(lonu,1)';
    latv = flipdim(latu,1)';
    lonu = lonu_tmp;
    latu = latu_tmp;
   end

   [lon2,lat2] = rot_sphere(lon1,lat1,rot);
   [lonu,latu] = rot_sphere(lonu,latu,rot);
   [lonv,latv] = rot_sphere(lonv,latv,rot);
   [lone,late] = rot_sphere(lone,late,rot);

   [lon3,lat3] = tra_sphere(lon2,lat2,tra_lat);
   [lonu,latu] = tra_sphere(lonu,latu,tra_lat);
   [lonv,latv] = tra_sphere(lonv,latv,tra_lat);
   [lone,late] = tra_sphere(lone,late,tra_lat);

   lon4 = lon3 + tra_lon*pi/180;
   lonu = lonu + tra_lon*pi/180;
   lonv = lonv + tra_lon*pi/180;
   lone = lone + tra_lon*pi/180;
   lon4(lon4<-pi) = lon4(lon4<-pi) + 2*pi;
   lonu(lonu<-pi) = lonu(lonu<-pi) + 2*pi;
   lonv(lonv<-pi) = lonv(lonv<-pi) + 2*pi;
   lone(lone<-pi) = lone(lone<-pi) + 2*pi;
   lat4 = lat3;

 %% Compute pn and pm
   %% pm = 1/dx
   pmu = gc_dist(lonu(:,1:end-1),latu(:,1:end-1),lonu(:,2:end),latu(:,2:end));
   pm = 0.*lon4;
   pm(:,2:end-1) = pmu;
   pm(:,1  ) = pm(:,2    );
   pm(:,end) = pm(:,end-1); 
   pm = 1./pm;

   %% pn = 1/dy
   pnv = gc_dist(lonv(1:end-1,:),latv(1:end-1,:),lonv(2:end,:),latv(2:end,:));
   pn = 0.*lon4;
   pn(2:end-1,:) = pnv;
   pn(1  ,:) = pn(2    ,:);
   pn(end,:) = pn(end-1,:);
   pn = 1./pn;


 %% Compute angles of local grid positive x-axis relative to east
   dellat = latu(:,2:end)-latu(:,1:end-1);
   dellon = lonu(:,2:end)-lonu(:,1:end-1);
   dellon(dellon> pi) = dellon(dellon>pi) - 2*pi;
   dellon(dellon<-pi) = dellon(dellon<-pi)+ 2*pi;
   dellon = dellon.*cos(0.5*(latu(:,2:end)+latu(:,1:end-1)));

   ang = lon4;
   ang_s = atan(dellat./(dellon+1e-16));
   ang_s(dellon<0 & dellat< 0) = ang_s(dellon<0 & dellat< 0) - pi;
   ang_s(dellon<0 & dellat>=0) = ang_s(dellon<0 & dellat>=0) + pi;
   ang_s(ang_s > pi) = ang_s(ang_s > pi) - pi;
   ang_s(ang_s <-pi) = ang_s(ang_s <-pi) + pi;

   ang(:,2:end-1) = ang_s;
   ang(:,1)   = ang(:,2);
   ang(:,end) = ang(:,end-1);

   lon4(lon4<0) = lon4(lon4<0) + 2*pi;
   lone(lone<0) = lone(lone<0) + 2*pi;

   return
