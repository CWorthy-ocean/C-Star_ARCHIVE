function [u,v,temp,salt,ssh,zi,loni,lati] = get_glorys_data(soda_mon_tr,soda_mon_vel,soda_mon_ssh,month,lon,lat);

% Processes and returns temp, salt and ssh from GLORYS

  disp('get_glorys')
% Get 3d temperature, salinity and mask from GLORYS
  soda_mon_data = soda_mon_tr;

  lon0 = min(lon(:));
  lon1 = max(lon(:));
  lat0 = min(lat(:));
  lat1 = max(lat(:));

  loni = ncread(soda_mon_data,'longitude');
  loni(loni<=0) = loni(loni<=0)+360;
  lati = ncread(soda_mon_data,'latitude');

  %%%%%%%%%%%
  %Extend domain if needed
  ext_east = 0;
  ext_west = 0;
  if lon0<min(loni)
    ext_west = 1    %% extend data on the west
    i0 = find(loni<lon0+360,1,'last');
  else
    i0 = find(loni<lon0,1,'last');
  end
  if lon1> max(loni)
    ext_east = 1    %% extend data on the east
    i1 = find(loni>lon1-360,1,'first')+1;
  else
    i1 = find(loni>lon1,1,'first')+1;
  end

  if ext_east & ext_west
    error 'extending twice'
  end

  if ext_east
    loni = [loni(i0:end)' loni(1:i1)'+360]';
  elseif ext_west
    loni = [loni(i0:end)'-360 loni(1:i1)']';
  else
    loni = loni(i0:i1);
  end

  %%%%%%%%%%%
  j0 = find(lati<lat0,1,'last')-1;
  j1 = find(lati>lat1,1,'first')+1;
  tnx = i1-i0+1;
  tny = j1-j0+1;

  lati = lati(j0:j1);
  [loni,lati] = meshgrid(loni,lati);

% figure(1);plot(lon,lat,'.k')
% figure(2);plot(loni,lati,'.k')
% error

%%%%%%%%%%%%%%%%%%%%%%%

  lev = ncread(soda_mon_data,'depth')';
  nz = length(lev);
% [month j0 j1 i0 i1]
  [ny,nx] = size(loni);
  zi = zeros(nz,ny,nx);
  for k = 1:nz
    zi(k,:,:) = -lev(k);
  end
  zi=flipud(zi);

%%%%%%%%%%%%%%%%%%%%%%%

   if ext_east|ext_west

    temp = pagetranspose([pagetranspose(ncread(soda_mon_tr,'thetao'  ,[i0 j0 1 1],[inf tny inf 1]))...
                          pagetranspose(ncread(soda_mon_tr,'thetao'  ,[1 j0 1 1],[i1 tny inf 1])) ]);
    salt = pagetranspose([pagetranspose(ncread(soda_mon_tr,'so'  ,[i0 j0 1 1],[inf tny inf 1]))...
                          pagetranspose(ncread(soda_mon_tr,'so'  ,[1 j0 1 1],[i1 tny inf 1])) ]);
    u = pagetranspose([pagetranspose(ncread(soda_mon_vel,'uo'  ,[i0 j0 1 1],[inf tny inf 1]))...
                       pagetranspose(ncread(soda_mon_vel,'uo'  ,[1 j0 1 1],[i1 tny inf 1])) ]);
    v = pagetranspose([pagetranspose(ncread(soda_mon_vel,'vo'  ,[i0 j0 1 1],[inf tny inf 1]))...
                       pagetranspose(ncread(soda_mon_vel,'vo'  ,[1 j0 1 1],[i1 tny inf 1])) ]);
    ssh = pagetranspose([pagetranspose(ncread(soda_mon_ssh,'zos'  ,[i0 j0 1],[inf tny 1]))...
                         pagetranspose(ncread(soda_mon_ssh,'zos'  ,[1 j0 1],[i1 tny 1])) ]);

   else

    temp = squeeze(ncread(soda_mon_tr,'thetao', [i0,j0,1,1],[tnx,tny,inf,1]));
    salt = squeeze(ncread(soda_mon_tr,'so', [i0,j0,1,1],[tnx,tny,inf,1]));
    u = squeeze(ncread(soda_mon_vel,'uo', [i0,j0,1,1],[tnx,tny,inf,1]));
    v = squeeze(ncread(soda_mon_vel,'vo', [i0,j0,1,1],[tnx,tny,inf,1]));
    ssh = squeeze(ncread(soda_mon_ssh,'zos', [i0,j0,1],[tnx,tny,1]));

   end

%%%%%%%%%%%%%%%%%%%%%%%

  disp('WARNING : GLORYS fix : copying deepest value below bottom')
  for i=1:i1-i0+1
    for j=1:j1-j0+1
      indnan = min(find(isnan(squeeze(temp(i,j,:))))) ;
      if (indnan~=1)
        temp(i,j,indnan:nz) = temp(i,j,indnan-1);
        salt(i,j,indnan:nz) = salt(i,j,indnan-1);
        u   (i,j,indnan:nz) = u   (i,j,indnan-1);
        v   (i,j,indnan:nz) = v   (i,j,indnan-1);
      end
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%

  temp(isnan(temp))=0;
  salt(isnan(salt))=0;
  u(isnan(u))=0;
  v(isnan(v))=0;
  ssh(isnan(ssh))=0;

%%%%%%%%%%%%%%%%%%%%%%%

  temp=flipud(permute(temp,[3,2,1]));
  salt=flipud(permute(salt,[3,2,1]));
  u=flipud(permute(u,[3,2,1]));
  v=flipud(permute(v,[3,2,1]));
  ssh=permute(ssh,[2,1]);

%%%%%%%%%%%%%%%%%%%%%%%

end

