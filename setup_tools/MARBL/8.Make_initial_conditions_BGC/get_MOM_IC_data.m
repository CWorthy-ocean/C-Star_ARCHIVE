function [var,zi,loni,lati] = get_MOM_IC_data(infile,variable,lon,lat) ;

momgridfile='/glade/work/mlevy/cesm_inputdata/ecosys_jan_IC_omip_MOM_tx0.66v1_c230310.nc'

%%%% Read lon/lat %%%%%%
%  lon1d = ncread(momgridfile,'NLON')';

%% WARNING: The longitude in the MOM6 global solution runs from around
% -290 E to 70 E. This is a simple fix to make this script for most of the world,
% but it will not work for domains between 0 and 70 E.

%  lon1d = lon1d + 360;
%  lat1d = ncread(momgridfile,'NLAT')';


  loni=ncread(momgridfile,'LON')';
  loni=loni+360;
  lati=ncread(momgridfile,'LAT')';

  i0_tmp = zeros(1,size(lati,1));
  i1_tmp = zeros(1,size(lati,1));
  for i = 1:size(lati,1)
    lon1d = squeeze(loni(i,:));
    i0_tmp(i) = find(lon1d<min(min(lon)),1,'last');
    i1_tmp(i) = find(lon1d>max(max(lon)),1,'first');
  end
  i0 = min(i0_tmp);
  i1 = max(i1_tmp);


  j0_tmp = zeros(1,i1-i0+1);
  j1_tmp = zeros(1,i1-i0+1);
  for j = i0:i1
    lat1d = squeeze(lati(:,j));
    j0_tmp(j-i0+1) = find(lat1d<min(min(lat)),1,'last');
    j1_tmp(j-i0+1) = find(lat1d>max(max(lat)),1,'first');

  end
  j0 = min(j0_tmp);
  j1 = max(j1_tmp);

  loni = loni(j0:j1,i0:i1);
  lati = lati(j0:j1,i0:i1);

%{

%momgridfile='/glade/work/mlevy/cesm_inputdata/MOM_IC.nc';

loni=ncread(momgridfile,'lonh');
loni=loni+360;
lati=ncread(momgridfile,'lath');

i0 = 300;
i1 = 450;
j0 = 300;
j1 = 450;
[loni,lati] = meshgrid(loni,lati);

loni = loni(i0:i1,j0:j1);
lati = lati(i0:i1,j0:j1);
%}

%  min(min(lon))
%  max(max(lon))
%  min(min(lat))
%  max(max(lat))
%  'yay'
%  min(min(loni))
%  max(max(loni))
%  min(min(lati))
%  max(max(lati))

  lev       = ncread(infile,'Layer')';
  nz = length(lev);

  [ny,nx] = size(loni);
  zi = zeros(nz,ny,nx);
  for k = 1:nz
    zi(k,:,:) = -lev(k);
  end

  var  = ncread(infile,variable);
%  disp('WARNING : GLORYS fix : copying deepest value below bottom')
%  for j=j0:j1
%    for i=i0:i1
%      indnan = min(find(isnan(squeeze(var(j,i,:))))) ;
%      if (indnan~=1)
%        var(i,j,indnan:nz) = var(i,j,indnan-1);
%      end
%    end
%  end
  var(isnan(var))=0;

  %var = var(j0:j1,i0:i1,:);
 % var = permute(squeeze(var),[3 2 1]);
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



