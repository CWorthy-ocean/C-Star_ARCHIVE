function s2r_hv_inibgc(BGC_INI,grdname,ininame,chdscd);

% Get S-coordinate params for child grid
  theta_s = chdscd.theta_s
  theta_b = chdscd.theta_b
  hc      = chdscd.hc
  Nc      = chdscd.N;
  scoord_c  = chdscd.scoord;

  h    = ncread(grdname,'h'       )';
  mask = ncread(grdname,'mask_rho')';
  lon  = ncread(grdname,'lon_rho' )';
  lat  = ncread(grdname,'lat_rho' )';
  lon(lon<0) = lon(lon<0) + 360;

  [Mc,Lc] = size(mask);
  maskc3d = zeros(Nc,Mc,Lc);
  for k = 1:Nc
    maskc3d(k,:,:) = mask;
  end

  % Z-coordinate (3D) on child grid
  zr = zlevs4(h, h*0, theta_s, theta_b, hc, Nc, 'r', scoord_c);
  zw = zlevs4(h, h*0, theta_s, theta_b, hc, Nc, 'w', scoord_c);

  for trc=1:length(BGC_INI.bgc_tracer)

    [var,zi,loni,lati] = get_MOM_data(BGC_INI.source{trc},BGC_INI.bgc_tracer{trc},lon,lat);
    input_var_name{trc} = BGC_INI.bgc_tracer{trc};

    tic
    disp('Computing interpolation coefficients');
    dummy_mask = squeeze(var(1,:,:))*0+1 ;
    [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
    A = get_hv_coef( zi, zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
    toc

    [nz,ny,nx] = size(var);
    var = reshape(A*reshape(var,nz*ny*nx,1),Nc,Mc,Lc);

    disp(' Writing ini file')
    ncwrite(ininame,input_var_name{trc},permute(var ,[3 2 1]),[1 1 1 1])
  end

return
