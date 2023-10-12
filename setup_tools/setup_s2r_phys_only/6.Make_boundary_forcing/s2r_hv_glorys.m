function s2r_hv_glorys(soda_mon_tr,soda_mon_vel,soda_mon_ssh,grdname,bryname,month,chdscd,obcflag);
%--------------------------------------------------------------
%  Produce a ROMS boundary file from SODA 2.0.4 data
%
%  Inspired by Roms_tools (IRD).
%  Thanks to Pierrick, Patrick (IRD), Francois (UCLA), Yusuke (UCLA)
%  Jeroen Molemaker (UCLA); nmolem@ucla.edu
%--------------------------------------------------------------

% Get S-coordinate params for child grid
  theta_b_c = chdscd.theta_b;
  theta_s_c = chdscd.theta_s;
  hc_c      = chdscd.hc;
  N_c       = chdscd.N;
  scoord_c  = chdscd.scoord;

%   nc = netcdf(grdname, 'nowrite');
%   [mpc npc] = size(nc{'h'}(:));
%   close(nc)
  [npc mpc] = size(ncread(grdname,'h'));

%   nd   = netcdf(bryname, 'write');

% Set bry_time
%   nd{'bry_time'}(month) = 15 + (month-1)*30;

 for bnd = 1:4
  disp('-------------------------------------------------------------')
  if ~obcflag(bnd)
    disp('Closed boundary')
    continue
  end
  if bnd==1
   disp('South boundary')
   i0 = 1;
   i1 = npc;
   j0 = 1;
   j1 = 2;
   fcoef = 'r2r_Glorys12_coefs_south.mat';
  end
  if bnd==2
   disp('East boundary')
   i0 = npc-1;
   i1 = npc;
   j0 = 1;
   j1 = mpc;
   fcoef = 'r2r_Glorys12_coefs_east.mat';
  end
  if bnd==3
   disp('North boundary')
   i0 = 1;
   i1 = npc;
   j0 = mpc-1;
   j1 = mpc;
   fcoef = 'r2r_Glorys12_coefs_north.mat';
  end
  if bnd==4
   disp('West boundary')
   i0 = 1;
   i1 = 2;
   j0 = 1;
   j1 = mpc;
   fcoef = 'r2r_Glorys12_coefs_west.mat';
  end

% Get topography data from childgrid

  hc    = ncread(grdname,'h'       ,[i0 j0],[i1-i0+1 j1-j0+1])';
  mask = ncread(grdname,'mask_rho',[i0 j0],[i1-i0+1 j1-j0+1])';
  angc = ncread(grdname,'angle'   ,[i0 j0],[i1-i0+1 j1-j0+1])';
  lon  = ncread(grdname,'lon_rho' ,[i0 j0],[i1-i0+1 j1-j0+1])';
  lat  = ncread(grdname,'lat_rho' ,[i0 j0],[i1-i0+1 j1-j0+1])';
  %lon(lon<0) = lon(lon<0) + 360;
  cosc  = cos(angc);         sinc  = sin(angc);

%plot(loni,lati,'.k')
%plot(lon,lat,'.r');

  [Mc,Lc] = size(mask);
  maskc3d = zeros(N_c,Mc,Lc);
  for k = 1:N_c
   maskc3d(k,:,:) = mask;
  end
  umask = maskc3d(:,:,2:end).*maskc3d(:,:,1:end-1);
  vmask = maskc3d(:,2:end,:).*maskc3d(:,1:end-1,:);

  % Z-coordinate (3D) on child grid
  zr = zlevs4(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'r', scoord_c);
  zw = zlevs4(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'w', scoord_c);
  [Nc Mc Lc] = size(zr);

%get_soda_data
%  [u,v,temp,salt,ssh,zi,loni,lati] = get_soda_data(soda_mon_tr,soda_mon_vel,soda_mon_ssh,1,lon,lat);
%get_glorys_data
  %[u,v,temp,salt,ssh,zi,loni,lati] = get_glorys_data(soda_mon_tr,soda_mon_vel,soda_mon_ssh,1,lon,lat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %if exist(fcoef,'file')
  %  disp(' ')
  %  disp('Reading interpolation coefficients from file');
  %  load(fcoef)
  %else
  soda_mon_data = soda_mon_tr{1};

  lon0 = min(lon(:));
  lon1 = max(lon(:));
  lat0 = min(lat(:));
  lat1 = max(lat(:));

  loni = ncread(soda_mon_data,'longitude');
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

  lev = ncread(soda_mon_tr{1},'depth')';
  nz = length(lev);
% [month j0 j1 i0 i1]
  [ny,nx] = size(loni);
  zi = zeros(nz,ny,nx);
  for k = 1:nz
    zi(k,:,:) = -lev(k);
  end
  zi=flipud(zi);

   if ext_east|ext_west
    ssh = pagetranspose([pagetranspose(ncread(soda_mon_ssh{1},'zos'  ,[i0 j0 1],[inf tny 1]))...
                         pagetranspose(ncread(soda_mon_ssh{1},'zos'  ,[1 j0 1],[i1 tny 1])) ]);
   else
    ssh = squeeze(ncread(soda_mon_ssh{1},'zos', [i0,j0,1],[tnx,tny,1]));
   end
    ssh(isnan(ssh))=0;
    ssh=permute(ssh,[2,1]);
    dummy_mask = ssh; dummy_mask(dummy_mask~=0)=1 ;




    tic
    disp(' ')
    disp('Computing interpolation coefficients');
    [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
%   plot(loni,lati,'.k');
%   hold on;plot(lon,lat,'.r');hold off
    A = get_hv_coef(zi, zr, coef2d, elem2d, loni, lati, lon, lat);
    save(fcoef,'elem2d','coef2d','nnel','A')
    toc
  %end

  for days = 1:length(soda_mon_tr)
     disp(['%%----> Working on days ' num2str(days) '/' num2str(length(soda_mon_tr)) ' <----%%'])

    [u,v,temp,salt,ssh,zi,loni,lati] = get_glorys_data(soda_mon_tr{days},soda_mon_vel{days},soda_mon_ssh{days},1,lon,lat);

    dummy_mask = ssh ; dummy_mask(dummy_mask~=0)=1 ;
  % fillmask
    temp=fillmask(temp, 1, dummy_mask, nnel);
    salt=fillmask(salt, 1, dummy_mask, nnel);

%   Prepare for estimating barotropic velocity
    dz  = zw(2:end,:,:)-zw(1:end-1,:,:);
    dzu = 0.5*(dz(:,:,1:end-1)+dz(:,:,2:end));
    dzv = 0.5*(dz(:,1:end-1,:)+dz(:,2:end,:));

%   Process scalar 3D variables
    for vint = 1:2 % Loop on the tracers
      if (vint==1)
        svar='temp';
	var = temp;
     	[nz,ny,nx] = size(temp);
      elseif (vint==2)
        svar='salt';
	var = salt;
      end
      var = reshape(A*reshape(var,nz*ny*nx,1),Nc,Mc,Lc);
      if (bnd ==1 )
        ncwrite(bryname,[svar '_south'],permute(var(:, 1, :),[3 1 2]),[1 1 days])
      end
      if (bnd ==2 )
        ncwrite(bryname,[svar '_east'],permute(var(:, :, end),[2 1 3]),[1 1 days])
      end
      if (bnd == 3 )
        ncwrite(bryname,[svar '_north'],permute(var(:, end, :),[3 1 2]),[1 1 days])
      end
      if (bnd == 4 )
        ncwrite(bryname,[svar '_west'],permute(var(:, :, 1),[2 1 3]),[1 1 days])
      end
    end  % End loop on vint

    % 3d interpolation of u_r and v_r to roms grid

    ud = reshape(A*reshape(u, nz*ny*nx,1), Nc,Mc,Lc);
    vd = reshape(A*reshape(v, nz*ny*nx,1), Nc,Mc,Lc);

    % Rotate to child orientation
    us = zeros(Nc,Mc,Lc);
    vs = zeros(Nc,Mc,Lc);
    for k=1:Nc
      us(k,:,:) = squeeze(ud(k,:,:)).*cosc + squeeze(vd(k,:,:)).*sinc;
      vs(k,:,:) = squeeze(vd(k,:,:)).*cosc - squeeze(ud(k,:,:)).*sinc;
    end
    u = 0.5*(us(:,:,1:Lc-1) + us(:,:,2:Lc));  %% back to staggered u points
    v = 0.5*(vs(:,1:Mc-1,:) + vs(:,2:Mc,:));  %% back to staggered v points

    u = u.*umask;
    v = v.*vmask;

    % Get barotropic velocity
    if sum(sum(sum(isnan(u)))) > 0
      error('nans in u velocity!')
    end
    if sum(sum(sum(isnan(v)))) > 0
      error('nans in v velocity!')
    end

%   disp('--- velocities');
%   disp(' Check this ')
    hu   = sum(dzu.*u); hv   = sum(dzv.*v);
    D_u  = sum(dzu);    D_v  = sum(dzv);
    [dum Mu Lu] = size(hu);
    [dum Mv Lv] = size(hv);
    ubar = reshape(hu./D_u,Mu, Lu);
    vbar = reshape(hv./D_v,Mv, Lv);


    % Sea surface height on ROMS grid
    zetac = sum(coef2d .* ssh(elem2d), 3);
    zetac = zetac.*mask;

%   disp(['------ saving velocity fields for record ' ...
%                       int2str(month) ' to bry file'])
    if (bnd ==1)
      ncwrite(bryname,'ubar_south',permute( ubar(1, :),[2 1]),[1 days])
      ncwrite(bryname,'vbar_south',permute( vbar(1, :),[2 1]),[1 days])
      ncwrite(bryname,'zeta_south',permute(zetac(1, :),[2 1]),[1 days])
      ncwrite(bryname,'u_south',permute(u(:, 1, :),[3 1 2]),[1 1 days])
      ncwrite(bryname,'v_south',permute(v(:, 1, :),[3 1 2]),[1 1 days])
    end
    if (bnd == 2)
      ncwrite(bryname,'ubar_east',         ubar(:,Lu)       ,[1 days])
      ncwrite(bryname,'vbar_east',         vbar(:,Lv)       ,[1 days])
      ncwrite(bryname,'zeta_east',        zetac(:,Lc)       ,[1 days])
      ncwrite(bryname,'u_east',permute(u(:,:,Lu),[2 1 3]),[1 1 days])
      ncwrite(bryname,'v_east',permute(v(:,:,Lv),[2 1 3]),[1 1 days])
    end
    if (bnd == 3)
      ncwrite(bryname,'ubar_north',permute( ubar(Mu, :),[2 1]),[1 days])
      ncwrite(bryname,'vbar_north',permute( vbar(Mv, :),[2 1]),[1 days])
      ncwrite(bryname,'zeta_north',permute(zetac(Mc, :),[2 1]),[1 days])
      ncwrite(bryname,'u_north',permute(u(:,Mu, :),[3 1 2]),[1 1 days])
      ncwrite(bryname,'v_north',permute(v(:,Mv, :),[3 1 2]),[1 1 days])
    end
    if (bnd == 4)
      ncwrite(bryname,'ubar_west',         ubar(:, 1)       ,[1 days])
      ncwrite(bryname,'vbar_west',         vbar(:, 1)       ,[1 days])
      ncwrite(bryname,'zeta_west',        zetac(:, 1)       ,[1 days])
      ncwrite(bryname,'u_west',permute(u(:, :, 1),[2 1 3]),[1 1 days])
      ncwrite(bryname,'v_west',permute(v(:, :, 1),[2 1 3]),[1 1 days])
    end

    % Close child bryfile

  end %days

 end    % End loop bnd

return























