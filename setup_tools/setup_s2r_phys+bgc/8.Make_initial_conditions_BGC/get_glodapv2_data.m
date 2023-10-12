function [var,zi,loni,lati] = get_glodapv2_data(glodapv2_file,variable,month,lon,lat,BRYtime,days);

% Processes and returns data from woa18 monthly climatology that is mask filled
% and extended from -180 to 360 degrees.

%%%% FIND the file to be read %%%%%%
  infile = glodapv2_file ;
  disp(['get_glodapv2 : reading file : ' infile])  
  
%%%% read lon/lat %%%%%%
  loni = ncread(infile,'lon')'; 
  disp('WARNING : special fix on longitude')
%   loni = [loni(181:end) loni(1:180)+360];
  lati = ncread(infile,'lat')';
  i0 = find(loni<min(min(lon)),1,'last');
  i1 = find(loni>max(max(lon)),1,'first');
  j0 = find(lati<min(min(lat)),1,'last');
  j1 = find(lati>max(max(lat)),1,'first');
  [loni,lati] = meshgrid(loni,lati);
% figure(1);plot(lon,lat,'.k')
% figure(2);plot(loni,lati,'.k')
% error
  loni = loni(j0:j1,i0:i1);
  lati = lati(j0:j1,i0:i1);
  
  lev = ncread(infile,'Depth')';
  nz = length(lev);
% [month j0 j1 i0 i1]
  [ny,nx] = size(loni);
  zi = zeros(nz,ny,nx);
  for k = 1:nz
    zi(k,:,:) = -lev(k);
  end
  
  if strcmp(variable,'DIC_glodap')==1
      name='TCO2'; factor = 1.0 ;
  elseif strcmp(variable,'Alk_glodap')==1
      name='TAlk'; factor = 1.0 ;
  else
      name=variable; factor = 1.0 ;
  end
  
  var = permute(squeeze(ncread(infile,name)),[3 2 1]);
  var = var(:,j0:j1,i0:i1)*factor ;
 
  disp('WARNING : fix for woa : bottom')
  for i=1:i1-i0+1
  for j=1:j1-j0+1
      indnan = min(find(isnan(squeeze(var(:,j,i))))) ;
      if (indnan~=1)
      var(indnan:nz,j,i) = var(indnan-1,j,i);
      end
  end
  end
  var(isnan(var))=0;  
  
  
disp('WARNING : return vertical in the ascending order') 
for k=1:nz
    var_new (k,:,:) = var (nz-k+1,:,:) ;
    zi_new  (k,:,:) = zi  (nz-k+1,:,:) ;
end
var = var_new ;
zi  = zi_new  ;

end

% mypcolor(loni,lati,squeeze(u(end,:,:)));set(gca,'ydir','normal');colorbar
% mypcolor(loni,lati,ssh);set(gca,'ydir','normal');colorbar
% load /batavia/nmolem/OBSERV/COAST/coast_atl
% hold on;plot(clon,clat,'k');hold off
% error('testing get soda data')


   
