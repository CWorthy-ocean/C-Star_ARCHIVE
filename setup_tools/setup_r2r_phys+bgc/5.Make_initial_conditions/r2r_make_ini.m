function r2r_make_ini(par_grd,par_data,chd_grd, chd_data,   ...
                      chdscd,   parscd, scoord_switch_p, scoord_switch_c );
%--------------------------------------------------------------
%
%  Make a roms 3d file for use as an initial file using data
%  from a roms parent grid.
%
%  input
%  =====
%  par_grd:         parent grid .nc file name
%  par_data:        parent history .nc file name
%  chd_grd:         child grid .nc file name
%  chd_data:        childt ini .nc file name
%  chdscd:          child grid s-coordinate parameters (object)
%  parscd:          parent grid s-coordinate parameters (object)
%  scoord_switch_c: child 'old' or 'new' s-coordinate
%  scoord_switch_p: parent 'old' or 'new' s-coordinate
%
%  Heavily modified from French produce
%  Jeroen Molemaker (UCLA); nmolem@atmos.ucla.edu
%--------------------------------------------------------------
%
%
% Get S-coordinate params for child grid
  N_c       = chdscd.N;
  theta_b_c = chdscd.theta_b;
  theta_s_c = chdscd.theta_s;
  hc_c      = chdscd.hc;

% Get S-coordinate params for parent grid
  N_p       = parscd.N;
  theta_b_p = parscd.theta_b;
  theta_s_p = parscd.theta_s;
  hc_p      = parscd.hc;
  par_file  = parscd.file;

% Get S-coordinate params for parent data file
  tind = parscd.tind;

% Set correct time in ini file
% np = netcdf(par_data, 'nowrite');
% ni = netcdf(chd_data, 'write');
  par_data
  chd_data
%   ncread(par_data,'time',tind,1)/(24*3600)
   ptime = ncread(par_data,'ocean_time',tind,1);
% % ni{'ocean_time'}(:) = np{'ocean_time'}(tind) + 100;  %% Adding 100 seconds to fix rounding problems in bry_time.
   ncwrite(chd_data,'ocean_time',ptime);
%  ncwrite(chd_data,'ocean_time',0);

% Get full parent grid and do triangulation
  lonp  = double(ncread(par_grd,'lon_rho')');%/10000;
  latp  = double(ncread(par_grd,'lat_rho')');%/10000;

  [Mpp,Lpp] = size(latp)
  lonp(lonp<0) = lonp(lonp<0) + 360;

  display('going delaunay');
% tri_fullpar = delaunay(lonp,latp);
  tri_fullpar = delaunayTriangulation(lonp(:),latp(:));
%  tri_fullpar = DelaunayTri([reshape(lonp,Mpp*Lpp,1),reshape(latp,Mpp*Lpp,1)]);
  display('return delaunay');

% Get child grid and chunk size
  ndomx = 7;
  ndomy = 7;
  [Mp,Lp] = size(ncread(chd_grd,'h')');
  szx = floor(Lp/ndomx);
  szy = floor(Mp/ndomy);

  icmin = [0:ndomx-1]*szx;
  jcmin = [0:ndomy-1]*szy;
  icmax = [1:ndomx]*szx;
  jcmax = [1:ndomy]*szy;
  icmin(1) = 1;
  jcmin(1) = 1;
  icmax(end) = Lp;
  jcmax(end) = Mp;

% Do the interpolation for all child chunks
  for domx = 1:ndomx
   for domy = 1:ndomy
     [ domx domy]

      icb = icmin(domx);
      ice = icmax(domx);
      jcb = jcmin(domy);
      jce = jcmax(domy);

    % Get topography data from childgrid
%     hc1 = ncread(chd_grd,'h')'; hc1 = hc1(jcb:jce,icb:ice);
      hc    = ncread(chd_grd,'h'       ,[icb jcb],[ice-icb+1 jce-jcb+1])';
      maskc = ncread(chd_grd,'mask_rho',[icb jcb],[ice-icb+1 jce-jcb+1])';
      angc  = ncread(chd_grd,'angle'   ,[icb jcb],[ice-icb+1 jce-jcb+1])';
      lonc  = double(ncread(chd_grd,'lon_rho' ,[icb jcb],[ice-icb+1 jce-jcb+1])');
      latc  = double(ncread(chd_grd,'lat_rho' ,[icb jcb],[ice-icb+1 jce-jcb+1])');
      umask = maskc(:,1:end-1).*maskc(:,2:end);
      vmask = maskc(1:end-1,:).*maskc(2:end,:);
      cosc  = cos(angc);
      sinc  = sin(angc);
      lonc(lonc<0) = lonc(lonc<0) + 360;
%     figure
%     plot(lonp,latp,'.k');
%     hold on
%     plot(lonc(1,:),latc(1,:),'.r');
%     plot(lonc(:,end),latc(:,end),'.r');
%     plot(lonc(end,:),latc(end,:),'.r');
%     plot(lonc(:,1),latc(:,1),'.r');
%     hold off
%     error 'testing'

    % Compute minimal subgrid extracted from full parent grid
%    t = squeeze(tsearch(lonp,latp,tri_fullpar,lonc,latc));
%       [nyc,nxc] = size(lonc);
%       t   = squeeze(pointLocation(tri_fullpar,reshape(lonc,nxc*nyc,1),reshape(latc,nxc*nyc,1)));
%       sum(isnan(t))

 [Mc,Lc] = size(lonc)

  xp = lonp; yp = latp;
  xc = lonc; yc = latc;

  Xp    = [reshape(xp,Mpp*Lpp,1) reshape(yp,Mpp*Lpp,1) ];
  Xc    = [reshape(xc,Mc*Lc,1) reshape(yc,Mc*Lc,1) ];

  ID = pointLocation(tri_fullpar,Xc);   % ID



    % Deal with child points that are outside parent grid (those points should be masked!)
      if (length(t(~isfinite(t)))>0);
%       disp('Warning in new_bry_subgrid: outside point(s) detected.');
%       [lonc,latc] = fix_outside_child(lonc,latc,t);
%       t = squeeze(tsearch(lonp,latp,tri_fullpar,lonc,latc));
       disp('ERROR in new_bry_subgrid: outside point(s) detected. Aborting.');
       return;
%       t = squeeze(pointLocation(tri_fullpar,reshape(double(lonc),nxc*nyc,1),reshape(double(latc),nxc*nyc,1)));
      end;
      %index       = tri_fullpar(t,:);
      index       = tri_fullpar(ID,:);
      [idxj,idxi] = ind2sub([Mpp Lpp], index);

      imin = min(min(idxi));% imin = max(1,imin-1);
      imax = max(max(idxi));% imax = min(1,imin-1);
      jmin = min(min(idxj));
      jmax = max(max(idxj));

    % Get parent grid and squeeze minimal subgrid
      [imin imax jmin jmax]
      par_grd
      masks = ncread(par_grd,'mask_rho',[imin jmin],[imax-imin+1,jmax-jmin+1])';
      lons  = ncread(par_grd,'lon_rho' ,[imin jmin],[imax-imin+1,jmax-jmin+1])'; %lons = double(lons)/1e4;
      lats  = ncread(par_grd,'lat_rho' ,[imin jmin],[imax-imin+1,jmax-jmin+1])'; %lats = double(lats)/1e4;
      angs  = ncread(par_grd,'angle'   ,[imin jmin],[imax-imin+1,jmax-jmin+1])';
      hs    = ncread(par_grd,'h'       ,[imin jmin],[imax-imin+1,jmax-jmin+1])';
      lons(lons<0) = lons(lons<0) + 360;
      coss = cos(angs); sins = sin(angs);
      if sum(isnan(masks))>0
        disp('Setting NaNs in masks to zero')
        error
        masks(isnan(masks))=0;
        disp('You probably have land masking defined in cppdefs.h...')
      end
%
%       figure
%       plot(lons,lats,'.k')
%       hold on
%       plot(lonc,latc,'.r')
%       hold off
%       return
%     size(hs)
    % Z-coordinate (3D) on minimal subgrid and child grid
      zs = zlevs4(hs, hs*0, theta_s_p, theta_b_p, hc_p, N_p, 'r', scoord_switch_p);
%     zs = zlev_cf(hs, hs*0, par_file, 'r');
      zc = zlevs4(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'r', scoord_switch_c);
      zw = zlevs4(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'w', scoord_switch_c);

      [Np Mp Lp] = size(zs);
      [Nc Mc Lc] = size(zc);

      disp('Computing interpolation coefficients');
%     lonc(lonc<0) = lonc(lonc<0) + 360;
%     plot(lons,lats,'.k');
%     hold on
%     plot(lonc(1,1),latc(1,1),'.r');
%       plot(lonc(1,end),latc(1,end),'.r');
%       plot(lonc(end,1),latc(end,1),'.r');
%       plot(lonc(end,end),latc(end,end),'.r');
%     hold off
%     error 'testing'
      [elem2d,coef2d,nnel] = get_tri_coef(lons,lats,lonc,latc,masks);
      A = get_hv_coef(zs, zc, coef2d, elem2d, lons, lats, lonc, latc);

    % Open parent data file


      disp('--- zeta')
      zetas = ncread(par_data,'zeta'  ,[imin jmin tind],[imax-imin+1,jmax-jmin+1 1])';
      zetas = fillmask(zetas, 1, masks, nnel);
      zetac = sum(coef2d.*zetas(elem2d), 3);
      zetac     = zetac.*maskc;

      disp(['--- temp'])
      var = squeeze(ncread(par_data,'temp'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
%     size(A)
%     save 'testing' 'A' 'var'
      ini_temp = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- salt'])
      var = squeeze(ncread(par_data,'salt'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_salt = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

    % Read in staggered velocities
      disp('--- baroclinic velocity');
      ud = squeeze(ncread(par_data,'u'  ,[imin jmin 1 tind],[imax-imin,jmax-jmin+1 inf 1]));
      vd = squeeze(ncread(par_data,'v'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin inf 1]));
      ud = permute(ud,[3 2 1]);
      vd = permute(vd,[3 2 1]);

      ud(abs(ud)>10) = 0;
      vd(abs(vd)>10) = 0;

    % Average to rho points
      ur = u2rho_3d(ud);
      vr = v2rho_3d(vd);

      if 0 %% Only for NHMG
    % compute vertical flux on parent grid
      hse = ncread(par_grd,'h',[imin-1 jmin-1],[imax-imin+3,jmax-jmin+3])';
      pms = ncread(par_grd,'pm',[imin-1 jmin-1],[imax-imin+3,jmax-jmin+3])';
      pns = ncread(par_grd,'pn',[imin-1 jmin-1],[imax-imin+3,jmax-jmin+3])';
      zws = zlevs3(hse, hse*0, theta_s_p, theta_b_p, hc_p, N_p, 'w', scoord_switch_p);
      dzs = zws(2:end,:,:)-zws(1:end-1,:,:);
      dzus = 0.5*(dzs(:,2:end-1,1:end-1)+dzs(:,2:end-1,2:end));
      dzvs = 0.5*(dzs(:,1:end-1,2:end-1)+dzs(:,2:end,2:end-1));

      pnus = 0.5*(pns(2:end-1,1:end-1)+pns(2:end-1,2:end));
      pmvs = 0.5*(pms(1:end-1,2:end-1)+pms(2:end,2:end-1));

      Flxu = ud;
      Flxv = vd;

      size(Flxu)
      size(dzus)
      size(ud)
      size(pnus)

      for k = 1:Np
       Flxu(k,:,:) = squeeze(dzus(k,:,:).*ud(k,:,:))./pnus;
       Flxv(k,:,:) = squeeze(dzvs(k,:,:).*vd(k,:,:))./pmvs;
      end

      ws = zeros(Np+1,Mp,Lp);
      for k = 1:Np
	    ws(k+1,:,:) = squeeze(ws(k,:,:)) - pms(2:end-1,2:end-1).*pns(2:end-1,2:end-1).*squeeze( ...
                Flxu(k,:,2:end)-Flxu(k,:,1:end-1) +  ...
                Flxv(k,2:end,:)-Flxv(k,1:end-1,:) );
      end
      % average ws to rho points
       ws = 0.5*(ws(2:end,:,:)+ws(1:end-1,:,:));
       ws = double(ws);
       wd = reshape(A*reshape(ws, Np*Mp*Lp,1), Nc,Mc,Lc);
      end

    % Rotate to north
      us = zeros(Np,Mp,Lp);
      vs = zeros(Np,Mp,Lp);
      for k = 1:Np
        us(k,:,:) = squeeze(ur(k,:,:)).*coss - squeeze(vr(k,:,:)).*sins;
        vs(k,:,:) = squeeze(vr(k,:,:)).*coss + squeeze(ur(k,:,:)).*sins;
      end

      us = fillmask(us, 0, masks, nnel);
      vs = fillmask(vs, 0, masks, nnel);
      us = double(us);
      vs = double(vs);
      ud = reshape(A*reshape(us, Np*Mp*Lp,1), Nc,Mc,Lc);
      vd = reshape(A*reshape(vs, Np*Mp*Lp,1), Nc,Mc,Lc);

    % Rotate to child orientation
      us = zeros(Nc, Mc, Lc);
      vs = zeros(Nc, Mc, Lc);
      for k=1:Nc
        us(k,:,:) = squeeze(ud(k,:,:)).*cosc + squeeze(vd(k,:,:)).*sinc;
        vs(k,:,:) = squeeze(vd(k,:,:)).*cosc - squeeze(ud(k,:,:)).*sinc;
      end
    % Back to staggered locations
      u = 0.5*(us(:,:,1:Lc-1) + us(:,:,2:Lc));
      v = 0.5*(vs(:,1:Mc-1,:) + vs(:,2:Mc,:));
    % w back to staggered w points
      if 0
      w = zeros(Nc+1,Mc,Lc);
      w(2:end-1,:,:) = 0.5*(wd(2:end,:,:)+wd(1:end-1,:,:));
      w(end,:,:) = w(end-1,:,:);
      end

    % Get barotropic velocity
      disp('--- barotropic velocity');
      dz  = zw(2:end,:,:)-zw(1:end-1,:,:);
      dzu = 0.5*(dz(:,:,1:end-1)+dz(:,:,2:end));
      dzv = 0.5*(dz(:,1:end-1,:)+dz(:,2:end,:));
      hu   = sum(dzu.*u); hv   = sum(dzv.*v);
      D_u  = sum(dzu);    D_v  = sum(dzv);
      ubar = squeeze(hu./D_u);     vbar = squeeze(hv./D_v);
      ubar      = ubar.*umask;
      vbar      = vbar.*vmask;

    % Zero-ing out the mask
      for k = 1:Nc
        ini_temp(k,:,:) = squeeze(ini_temp(k,:,:)); %.*maskc;
        ini_salt(k,:,:) = squeeze(ini_salt(k,:,:)); %.*maskc;
        u(k,:,:)        = squeeze(u(k,:,:)).*umask;
        v(k,:,:)        = squeeze(v(k,:,:)).*vmask;
      end

%%%%%%%%%%%%%%%%

      disp(['--- PO4'])
      var = squeeze(ncread(par_data,'PO4'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_PO4 = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- NO3'])
      var = squeeze(ncread(par_data,'NO3'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_NO3 = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- SiO3'])
      var = squeeze(ncread(par_data,'SiO3'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_SiO3 = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- NH4'])
      var = squeeze(ncread(par_data,'NH4'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_NH4 = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- Fe'])
      var = squeeze(ncread(par_data,'Fe'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_Fe = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- O2'])
      var = squeeze(ncread(par_data,'O2'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_O2 = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DIC'])
      var = squeeze(ncread(par_data,'DIC'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DIC = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- Alk'])
      var = squeeze(ncread(par_data,'Alk'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_Alk = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DOC'])
      var = squeeze(ncread(par_data,'DOC'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DOC = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DON'])
      var = squeeze(ncread(par_data,'DON'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DON = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DOFE'])
      var = squeeze(ncread(par_data,'DOFE'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DOFE = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DOP'])
      var = squeeze(ncread(par_data,'DOP'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DOP = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DOPR'])
      var = squeeze(ncread(par_data,'DOPR'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DOPR = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DONR'])
      var = squeeze(ncread(par_data,'DONR'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DONR = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- ZOOC'])
      var = squeeze(ncread(par_data,'ZOOC'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_ZOOC = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- SPC'])
      var = squeeze(ncread(par_data,'SPC'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_SPC = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- SPCHL'])
      var = squeeze(ncread(par_data,'SPCHL'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_SPCHL = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- SPFE'])
      var = squeeze(ncread(par_data,'SPFE'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_SPFE = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- SPCACO3'])
      var = squeeze(ncread(par_data,'SPCACO3'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_SPCACO3 = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DIATC'])
      var = squeeze(ncread(par_data,'DIATC'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DIATC = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DIATCHL'])
      var = squeeze(ncread(par_data,'DIATCHL'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DIATCHL = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DIATFE'])
      var = squeeze(ncread(par_data,'DIATFE'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DIATFE = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DIATSI'])
      var = squeeze(ncread(par_data,'DIATSI'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DIATSI = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DIAZC'])
      var = squeeze(ncread(par_data,'DIAZC'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DIAZC = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DIAZCHL'])
      var = squeeze(ncread(par_data,'DIAZCHL'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DIAZCHL = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- DIAZFE'])
      var = squeeze(ncread(par_data,'DIAZFE'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_DIAZFE = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- NO2'])
      var = squeeze(ncread(par_data,'NO2'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_NO2 = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- N2'])
      var = squeeze(ncread(par_data,'N2'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_N2 = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);

      disp(['--- N2O'])
      var = squeeze(ncread(par_data,'N2O'  ,[imin jmin 1 tind],[imax-imin+1,jmax-jmin+1 inf 1]));
      var = permute(var,[3 2 1]);
      var = fillmask(var,1,masks,nnel);
      var = double(var);
      ini_N2O = reshape(A*reshape(var,Np*Mp*Lp,1),Nc,Mc,Lc);


      disp(' Writing ini file')
      ini_temp = permute(ini_temp,[3 2 1]);
      ncwrite(chd_data,'temp'  ,ini_temp,[icb jcb 1 1]);
      ini_salt = permute(ini_salt,[3 2 1]);
      ncwrite(chd_data,'salt'  ,ini_salt,[icb jcb 1 1]);
      u = permute(u,[3 2 1]);
      ncwrite(chd_data,'u'     ,u       ,[icb jcb 1 1]);
      v = permute(v,[3 2 1]);
      ncwrite(chd_data,'v'     ,v       ,[icb jcb 1 1]);
%       w = permute(w,[3 2 1]);
%       ncwrite(chd_data,'w'     ,w       ,[icb jcb 1 1]);
      ncwrite(chd_data,'zeta'  ,zetac'  ,[icb jcb 1]  );
%      ncwrite(chd_data,'zeta'  ,zetac'  ,[icb jcb 1]  );
      ncwrite(chd_data,'ubar'  ,ubar'   ,[icb jcb 1]  );
      ncwrite(chd_data,'vbar'  ,vbar'   ,[icb jcb 1]  );

      ini_PO4 = permute(ini_PO4,[3 2 1]);
      ncwrite(chd_data,'PO4'  ,ini_PO4,[icb jcb 1 1]);
      ini_NO3 = permute(ini_NO3,[3 2 1]);
      ncwrite(chd_data,'NO3'  ,ini_NO3,[icb jcb 1 1]);
      ini_SiO3 = permute(ini_SiO3,[3 2 1]);
      ncwrite(chd_data,'SiO3'  ,ini_SiO3,[icb jcb 1 1]);
      ini_NH4 = permute(ini_NH4,[3 2 1]);
      ncwrite(chd_data,'NH4'  ,ini_NH4,[icb jcb 1 1]);
      ini_Fe = permute(ini_Fe,[3 2 1]);
      ncwrite(chd_data,'Fe'  ,ini_Fe,[icb jcb 1 1]);
      ini_O2 = permute(ini_O2,[3 2 1]);
      ncwrite(chd_data,'O2'  ,ini_O2,[icb jcb 1 1]);
      ini_DIC = permute(ini_DIC,[3 2 1]);
      ncwrite(chd_data,'DIC'  ,ini_DIC,[icb jcb 1 1]);
      ini_Alk = permute(ini_Alk,[3 2 1]);
      ncwrite(chd_data,'Alk'  ,ini_Alk,[icb jcb 1 1]);
      ini_DOC = permute(ini_DOC,[3 2 1]);
      ncwrite(chd_data,'DOC'  ,ini_DOC,[icb jcb 1 1]);
      ini_DON = permute(ini_DON,[3 2 1]);
      ncwrite(chd_data,'DON'  ,ini_DON,[icb jcb 1 1]);
      ini_DOFE = permute(ini_DOFE,[3 2 1]);
      ncwrite(chd_data,'DOFE'  ,ini_DOFE,[icb jcb 1 1]);
      ini_DOP = permute(ini_DOP,[3 2 1]);
      ncwrite(chd_data,'DOP'  ,ini_DOP,[icb jcb 1 1]);
      ini_DOPR = permute(ini_DOPR,[3 2 1]);
      ncwrite(chd_data,'DOPR'  ,ini_DOPR,[icb jcb 1 1]);
      ini_DONR = permute(ini_DONR,[3 2 1]);
      ncwrite(chd_data,'DONR'  ,ini_DONR,[icb jcb 1 1]);
      ini_ZOOC = permute(ini_ZOOC,[3 2 1]);
      ncwrite(chd_data,'ZOOC'  ,ini_ZOOC,[icb jcb 1 1]);
      ini_SPC = permute(ini_SPC,[3 2 1]);
      ncwrite(chd_data,'SPC'  ,ini_SPC,[icb jcb 1 1]);
      ini_SPCHL = permute(ini_SPCHL,[3 2 1]);
      ncwrite(chd_data,'SPCHL'  ,ini_SPCHL,[icb jcb 1 1]);
      ini_SPFE = permute(ini_SPFE,[3 2 1]);
      ncwrite(chd_data,'SPFE'  ,ini_SPFE,[icb jcb 1 1]);
      ini_SPCACO3 = permute(ini_SPCACO3,[3 2 1]);
      ncwrite(chd_data,'SPCACO3'  ,ini_SPCACO3,[icb jcb 1 1]);
      ini_DIATC = permute(ini_DIATC,[3 2 1]);
      ncwrite(chd_data,'DIATC'  ,ini_DIATC,[icb jcb 1 1]);
      ini_DIATCHL = permute(ini_DIATCHL,[3 2 1]);
      ncwrite(chd_data,'DIATCHL'  ,ini_DIATCHL,[icb jcb 1 1]);
      ini_DIATFE = permute(ini_DIATFE,[3 2 1]);
      ncwrite(chd_data,'DIATFE'  ,ini_DIATFE,[icb jcb 1 1]);
      ini_DIATSI = permute(ini_DIATSI,[3 2 1]);
      ncwrite(chd_data,'DIATSI'  ,ini_DIATSI,[icb jcb 1 1]);
      ini_DIAZC = permute(ini_DIAZC,[3 2 1]);
      ncwrite(chd_data,'DIAZC'  ,ini_DIAZC,[icb jcb 1 1]);
      ini_DIAZCHL = permute(ini_DIAZCHL,[3 2 1]);
      ncwrite(chd_data,'DIAZCHL'  ,ini_DIAZCHL,[icb jcb 1 1]);
      ini_DIAZFE = permute(ini_DIAZFE,[3 2 1]);
      ncwrite(chd_data,'DIAZFE'  ,ini_DIAZFE,[icb jcb 1 1]);
      ini_NO2 = permute(ini_NO2,[3 2 1]);
      ncwrite(chd_data,'NO2'  ,ini_NO2,[icb jcb 1 1]);
      ini_N2 = permute(ini_N2,[3 2 1]);
      ncwrite(chd_data,'N2'  ,ini_N2,[icb jcb 1 1]);
      ini_N2O = permute(ini_N2O,[3 2 1]);
      ncwrite(chd_data,'N2O'  ,ini_N2O,[icb jcb 1 1]);
   end
  end
