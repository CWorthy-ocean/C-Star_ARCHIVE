function [var,zi,loni,lati] = get_MOM_data(infile,variable,lon,lat) ;

%%%% Read lon/lat %%%%%%
  loni = ncread(infile,'lonh')';

%% WARNING: The longitude in the MOM6 global solution runs from around
% -290 E to 70 E. This is a simple fix to make this script for most of the world,
% but it will not work for domains between 0 and 70 E.

  loni = loni + 360;
  lati = ncread(infile,'lath')';
  i0 = find(loni<min(min(lon)),1,'last');
  i1 = find(loni>max(max(lon)),1,'first');
  j0 = find(lati<min(min(lat)),1,'last');
  j1 = find(lati>max(max(lat)),1,'first');

  [loni,lati] = meshgrid(loni,lati);
  loni = loni(j0:j1,i0:i1);
  lati = lati(j0:j1,i0:i1);

  lev       = ncread(infile,'Layer')';
  nz = length(lev);

  [ny,nx] = size(loni);
  zi = zeros(nz,ny,nx);
  for k = 1:nz
    zi(k,:,:) = -lev(k);
  end

  var  = ncread(infile,variable);
  var = permute(squeeze(var),[3 2 1]);
  var = var(:,j0:j1,i0:i1);

%% Reverse the order of the vertical
for k=1:nz
  var_new (k,:,:) = var (nz-k+1,:,:) ;
  zi_new  (k,:,:) = zi  (nz-k+1,:,:) ;
end
var = var_new ;
zi  = zi_new  ;

end



