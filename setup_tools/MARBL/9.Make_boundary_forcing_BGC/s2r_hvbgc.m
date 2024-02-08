function s2r_hvbgc(BGC_INI,grdname,bryname,chdscd,obcflag,BRYtime);
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

  [npc mpc] = size(ncread(grdname,'h'));


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
  end
  if bnd==2
   disp('East boundary')
   i0 = npc-1;
   i1 = npc;
   j0 = 1;
   j1 = mpc;
  end
  if bnd==3
   disp('North boundary')
   i0 = 1;
   i1 = npc;
   j0 = mpc-1;
   j1 = mpc;
  end
  if bnd==4
   disp('West boundary')
   i0 = 1;
   i1 = 2;
   j0 = 1;
   j1 = mpc;
  end

% Get topography data from childgrid

  hc    = ncread(grdname,'h'       ,[i0 j0],[i1-i0+1 j1-j0+1])';
  mask = ncread(grdname,'mask_rho',[i0 j0],[i1-i0+1 j1-j0+1])';
  lon  = ncread(grdname,'lon_rho' ,[i0 j0],[i1-i0+1 j1-j0+1])';
  lat  = ncread(grdname,'lat_rho' ,[i0 j0],[i1-i0+1 j1-j0+1])';
%   lon = lon + 360;
  lon(lon<0) = lon(lon<0) + 360;

  [Mc,Lc] = size(mask);
  maskc3d = zeros(N_c,Mc,Lc);
  for k = 1:N_c
   maskc3d(k,:,:) = mask;
  end

  % Z-coordinate (3D) on child grid
  zr = zlevs4(hc, hc*0, theta_s_c, theta_b_c, hc_c, N_c, 'r', scoord_c);
  [Nc Mc Lc] = size(zr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for trc=1:length(BGC_INI.bgc_tracer)

  for t = 1:length(BRYtime.source)


    [var,zi,loni,lati] = get_MOM_BC_data(BRYtime.source{t},BGC_INI.bgc_tracer{trc},lon,lat);
    input_var_name{trc} = BGC_INI.bgc_tracer{trc};

    tic
    disp('Computing interpolation coefficients');
    dummy_mask = squeeze(var(1,:,:))*0+1 ;
    [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
    A = get_hv_coef( zi, zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
    toc

    [nz,ny,nx] = size(var);

    %   Process scalar 3D variables
    var = reshape(A*reshape(var,nz*ny*nx,1),Nc,Mc,Lc);

    var = inpaintn(var);

    if (bnd ==1 )
       ncwrite(bryname,[input_var_name{trc} '_south'],permute(var(:, 1, :),[3 1 2]),[1 1 t])
    end
    if (bnd ==2 )
       ncwrite(bryname,[input_var_name{trc} '_east'],permute(var(:, :, end),[2 1 3]),[1 1 t])
    end
    if (bnd == 3 )
       ncwrite(bryname,[input_var_name{trc} '_north'],permute(var(:, end, :),[3 1 2]),[1 1 t])
    end
    if (bnd == 4 )
       ncwrite(bryname,[input_var_name{trc} '_west'],permute(var(:, :, 1),[2 1 3]),[1 1 t])
    end

  end

end



end % which boundary

return
