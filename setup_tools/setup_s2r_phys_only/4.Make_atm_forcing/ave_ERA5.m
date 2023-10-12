
frc_dir = '/paracas/DATASETS/ERA5/';

 eralist = dir([frc_dir 'ERA5*']);
 nfiles = length(eralist);

  swr_av = zeros(1440,721,12);
  cnt = zeros(12,1);
  for i = 1:nfiles
    datname = [frc_dir eralist(i).name] 
    mstr = datname(end-4:end-3);
    month = str2num(mstr);

    swr = mean(ncread(datname,'ssr'),3);

    swr_av(:,:,month) = swr_av(:,:,month) + swr;
    cnt(month) = cnt(month) + 1;
  end

  for month = 1:12
    swr_av(:,:,month) = swr_av(:,:,month)/cnt(month);
  end

  remove('ERA_swr_clim.nc')
  nccreate('ERA_swr_clim.nc','swr_av','Dimensions',{'lon', 1440,'lat',721,'time',12})
  ncwrite('ERA_swr_clim.nc','swr_av',swr_av);

  lon = ncread(datname,'longitude');
  lat = ncread(datname,'latitude');

  nccreate('ERA_swr_clim.nc','longitude','Dimensions',{'lon', 1440})
  nccreate('ERA_swr_clim.nc','latitude','Dimensions',{'lat', 721})
  ncwrite('ERA_swr_clim.nc','latitude',lat);
  ncwrite('ERA_swr_clim.nc','longitude',lon);

