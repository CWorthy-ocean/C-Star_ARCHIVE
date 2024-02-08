function [var,zi,loni,lati] = get_MOM_data(infile,variable,lon,lat) ;

  momgridfile='/glade/work/mlevy/cesm_inputdata/ecosys_jan_IC_omip_MOM_tx0.66v1_c230310.nc';

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



  lev       = ncread(infile,'zl')';
  nz = length(lev);

  [ny,nx] = size(loni);
  zi = zeros(nz,ny,nx);
  for k = 1:nz
    zi(k,:,:) = -lev(k);
  end

  if strcmp(variable,'spCaCO3')
     var = ncread(infile,'coccoCaCO3');
  elseif strcmp(variable,'zooC')
     var = ncread(infile,'mesozooC') + ncread(infile,'microzooC');
  elseif strcmp(variable,'ALK_ALT_CO2')
     var = ncread(infile,'ALK');
  elseif strcmp(variable,'DIC_ALT_CO2')
     var = ncread(infile,'DIC');
  else
     var = ncread(infile,variable);
  end

  var(isnan(var))=0;

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



