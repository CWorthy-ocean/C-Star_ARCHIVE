 function dis = gc_dist(lon1,lat1,lon2,lat2)
%
%  Distance between 2 points along a great circle
%  lat and lon in radians!!
%  2008, Jeroen Molaker, UCLA

   dlat = lat2-lat1;
   dlon = lon2-lon1;
%  dlon(dlon>2*pi) = dlon(dlon>2*pi)-2*pi;
   dang = 2*asin( sqrt( sin(dlat/2).^2 + cos(lat2).*cos(lat1).*sin(dlon/2).^2 ) );  %% haversine function

   r_earth = 6371315.;

   dis = r_earth*dang;

