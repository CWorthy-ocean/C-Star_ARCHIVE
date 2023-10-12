function [u,v,temp,salt,ssh,zi,loni,lati] = get_soda_data(soda_mon_tr,soda_mon_vel,soda_mon_ssh,month,lon,lat);

% Processes and returns temp, salt and ssh from SODA monthly climatology that is mask filled
% and extended from -180 to 360 degrees.

  disp('get_soda')
% Get 3d temperature, salinity and mask from SODA climatology
  soda_mon_data = soda_mon_tr;
  loni = ncread(soda_mon_data,'xt_ocean')';
%   loni = ncread(soda_mon_data,'xt_ocean',1,1120)'; 
  disp('WARNING : special fix on longitude made for 5days soda to Pacifique')
%   loni = [loni(1:1120)] + 360 ; %%% SPECIAL FOR SODA 5days AVG
  loni = [loni(1121:1440) loni(1:1120)+ 360]  ; %%% SPECIAL FOR SODA 5days AVG
%   loni = loni + 360 ;
  %%%%%%%%% Ugly fix   %%%%%%%%%%%%%%%%%%%%%%%
%   loni(1) = 79.83 ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  lati = ncread(soda_mon_data,'yt_ocean')';
  i0 = find(loni<min(min(lon)),1,'last');
  if (isempty(i0))
      disp('WARNING 2: special fix on min longitude made for 5days soda to Pacifique')
      i0=1;
  end
  i1 = find(loni>max(max(lon)),1,'first')+1;
  j0 = find(lati<min(min(lat)),1,'last')-1;
  j1 = find(lati>max(max(lat)),1,'first')+1;
  [loni,lati] = meshgrid(loni,lati);
% figure(1);plot(lon,lat,'.k')
% figure(2);plot(loni,lati,'.k')
% error
  loni = loni(j0:j1,i0:i1);
  lati = lati(j0:j1,i0:i1);
  
  lev = ncread(soda_mon_data,'st_ocean')';
  nz = length(lev);
% [month j0 j1 i0 i1]
  [ny,nx] = size(loni);
  zi = zeros(nz,ny,nx);
  for k = 1:nz
    zi(k,:,:) = -lev(k);
  end
  

%   temp = permute(squeeze(ncread(soda_mon_data,'temp',[i0 j0 1 month],[i1-i0+1 j1-j0+1 Inf 1])),[3 2 1]);
%   salt = permute(squeeze(ncread(soda_mon_data,'salt',[i0 j0 1 month],[i1-i0+1 j1-j0+1 Inf 1])),[3 2 1]);
  temp = ncread(soda_mon_data,'temp') ;  
  temp = permute([permute(temp(1121:1440,:,:),[2 1 3]) permute(temp(1:1120,:,:),[2 1 3])],[2 1 3]) ;
  temp = temp(i0:i1,j0:j1,:) ; temp = permute(temp,[3 2 1]) ;
  salt = ncread(soda_mon_data,'salt') ; 
  salt = permute([permute(salt(1121:1440,:,:),[2 1 3]) permute(salt(1:1120,:,:),[2 1 3])],[2 1 3]) ;
  salt = salt(i0:i1,j0:j1,:) ; salt = permute(salt,[3 2 1]) ;
  
  soda_mon_data = soda_mon_vel;
%   u = permute(squeeze(ncread(soda_mon_data,'u',[i0 j0 1 month],[i1-i0+1 j1-j0+1 Inf 1])),[3 2 1]);
%   v = permute(squeeze(ncread(soda_mon_data,'v',[i0 j0 1 month],[i1-i0+1 j1-j0+1 Inf 1])),[3 2 1]);
  u = ncread(soda_mon_data,'u') ;  
  u = permute([permute(u(1121:1440,:,:),[2 1 3]) permute(u(1:1120,:,:),[2 1 3])],[2 1 3]) ;
  u = u(i0:i1,j0:j1,:) ; u = permute(u,[3 2 1]) ;  
  v = ncread(soda_mon_data,'v') ; 
  v = permute([permute(v(1121:1440,:,:),[2 1 3]) permute(v(1:1120,:,:),[2 1 3])],[2 1 3]) ;
  v = v(i0:i1,j0:j1,:) ; v = permute(v,[3 2 1]) ;  
  
  
  disp('WARNING : fix for 5days soda to Pacifique : bottom')
  for i=1:i1-i0+1
  for j=1:j1-j0+1
      indnan = min(find(isnan(squeeze(temp(:,j,i))))) ;
      if (indnan~=1)
      temp(indnan:nz,j,i) = temp(indnan-1,j,i);
      salt(indnan:nz,j,i) = salt(indnan-1,j,i);
      u   (indnan:nz,j,i) = u   (indnan-1,j,i);
      v   (indnan:nz,j,i) = v   (indnan-1,j,i);
      end
  end
  end
  temp(isnan(temp))=0;  
  salt(isnan(salt))=0;
  u(isnan(u))=0;
  v(isnan(v))=0;
  
  soda_mon_data = soda_mon_ssh ;  
%   ssh = permute(squeeze(ncread(soda_mon_data,'ssh',[i0 j0 month],[i1-i0+1 j1-j0+1 1])),[2 1]);
  ssh = ncread(soda_mon_data,'ssh') ; 
  ssh = [ssh(1121:1440,:)' ssh(1:1120,:)']'  ;
  ssh = ssh(i0:i1,j0:j1) ; ssh = permute(ssh,[2 1]) ;  
  ssh(isnan(ssh))=0;
  
disp('WARNING : return 5days soda vertical in the ascending order') 
for k=1:nz
    temp_new(k,:,:) = temp(nz-k+1,:,:) ;
    salt_new(k,:,:) = salt(nz-k+1,:,:) ;
    u_new   (k,:,:) = u   (nz-k+1,:,:) ;
    v_new   (k,:,:) = v   (nz-k+1,:,:) ; 
    zi_new  (k,:,:) = zi  (nz-k+1,:,:) ;
end
temp= temp_new ;
salt= salt_new ;
u   = u_new    ;
v   = v_new    ;
zi  = zi_new   ;

end

% mypcolor(loni,lati,squeeze(u(end,:,:)));set(gca,'ydir','normal');colorbar
% mypcolor(loni,lati,ssh);set(gca,'ydir','normal');colorbar
% load /batavia/nmolem/OBSERV/COAST/coast_atl
% hold on;plot(clon,clat,'k');hold off
% error('testing get soda data')


   
