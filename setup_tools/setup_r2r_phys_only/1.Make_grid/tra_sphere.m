function [lon2,lat2] = tra_sphere(lon1,lat1,tra)
%
%  Rotate sphere around its y-axis
%  Part of easy grid
%  (c) 2008, Jeroen Molemaker, UCLA

   
   [n,m] = size(lon1);
   tra = tra*pi/180;  %% translation in latitude direction

   % translate into x,y,z
   % conventions:  (lon,lat) = (0,0)  corresponds to (x,y,z) = ( 0,-r, 0)
   %               (lon,lat) = (0,90) corresponds to (x,y,z) = ( 0, 0, r)
   x1 = sin(lon1).*cos(lat1);
   y1 = cos(lon1).*cos(lat1);
   z1 = sin(lat1);

   % We will rotate these points around the small circle defined by 
   % the intersection of the sphere and the plane that
   % is orthogonal to the line through (lon,lat) (90,0) and (-90,0)

   % The rotation is in that plane around its intersection with
   % aforementioned line.

   % Since the plane is orthogonal to the x-axis (in my definition at least),
   % Rotations in the plane of the small circle maintain constant x and are around
   % (x,y,z) = (x1,0,0)

     rp1 = sqrt(y1.^2+z1.^2);

     ap1 = pi/2*ones(n,m);
     ap1(abs(y1)>1e-7) = atan(abs(z1(abs(y1)>1e-7)./y1(abs(y1)>1e-7)));
     ap1(y1<0) = pi - ap1(y1<0);
     ap1(z1<0) = -ap1(z1<0);

     ap2 = ap1 + tra;
     x2 = x1;
     y2 = rp1.*cos(ap2);
     z2 = rp1.*sin(ap2);

     %% transformation from (x,y,z) to (lat,lon)
     lon2 = pi/2*ones(n,m);
     lon2(abs(y2)>1e-7) = atan(abs(x2(abs(y2)>1e-7)./y2(abs(y2)>1e-7)));
     lon2(y2<0) = pi - lon2(y2<0);
     lon2(x2<0) = -lon2(x2<0);

     pr2 = sqrt(x2.^2+y2.^2);
     lat2 = pi/2*ones(n,m);
     lat2(abs(pr2)>1e-7) = atan(abs(z2(abs(pr2)>1e-7)./pr2(abs(pr2)>1e-7)));
     lat2(z2<0) = -lat2(z2<0);

   return
   
