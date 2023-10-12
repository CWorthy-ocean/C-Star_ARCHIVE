function s2r_hvbgc(BGC_INI,grdname,bryname,days,chdscd,obcflag,BRYtime);
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
%    if (BGC_INI.data_woa18==1)
%    fcoef_woa18 = 'r2r_woa18_coefs_south.mat';
%    end
  end
  if bnd==2
   disp('East boundary')
   i0 = npc-1;
   i1 = npc;
   j0 = 1;
   j1 = mpc;
%    if (BGC_INI.data_woa18==1)
%    fcoef_woa18 = 'r2r_Soda_coefs_east.mat';
%    end
  end
  if bnd==3
   disp('North boundary')
   i0 = 1;
   i1 = npc;
   j0 = mpc-1;
   j1 = mpc;
%    if (BGC_INI.data_woa18==1)
%    fcoef_woa18 = 'r2r_Soda_coefs_north.mat';
%    end
  end
  if bnd==4
   disp('West boundary')
   i0 = 1;
   i1 = 2;
   j0 = 1;
   j1 = mpc;
%    if (BGC_INI.data_woa18==1)
%    fcoef_woa18 = 'r2r_Soda_coefs_west.mat';
%    end
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
  if (BGC_INI.data_woa18==1)

      clear vars vars_new input_var_name
      for trc=1:length(BGC_INI.tracer_woa18)

            [var,zi,loni,lati] = get_woa18_data(BGC_INI.bgc_frcini{BGC_INI.tracer_woa18(trc)}, ...
                       BGC_INI.bgc_tracer{BGC_INI.tracer_woa18(trc)},BRYtime.month(days),lon,lat,BRYtime,days);
            input_var_name{trc} = BGC_INI.bgc_tracer{BGC_INI.tracer_woa18(trc)} ;

            dummy_mask = squeeze(var(end,:,:)) ; dummy_mask(dummy_mask~=0)=1 ;
%             if exist(fcoef_woa18,'file')
%                disp(' ')
%                disp('Reading interpolation coefficients from file');
%                load(fcoef_woa18)
%             else
               tic
               disp(' ')
               disp('Computing interpolation coefficients');
               [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati), ...
                          lon,lat,dummy_mask) ;
               A = get_hv_coef(zi, zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
%                save(fcoef_woa18,'elem2d','coef2d','nnel','A')
               toc
%             end

           % fillmask
           var=fillmask(var, 1, dummy_mask, nnel);

           [nz,ny,nx] = size(var);

           %   Process scalar 3D variables
           var = reshape(A*reshape(var,nz*ny*nx,1),Nc,Mc,Lc);
           if (bnd ==1 )
              ncwrite(bryname,[input_var_name{trc} '_south'],permute(var(:, 1, :),[3 1 2]),[1 1 days])
           end
           if (bnd ==2 )
              ncwrite(bryname,[input_var_name{trc} '_east'],permute(var(:, :, end),[2 1 3]),[1 1 days])
           end
           if (bnd == 3 )
              ncwrite(bryname,[input_var_name{trc} '_north'],permute(var(:, end, :),[3 1 2]),[1 1 days])
           end
           if (bnd == 4 )
              ncwrite(bryname,[input_var_name{trc} '_west'],permute(var(:, :, 1),[2 1 3]),[1 1 days])
           end

      end    % End loop trc

  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (BGC_INI.data_CCSM==1)

      clear vars vars_new input_var_name
      for trc=1:length(BGC_INI.tracer_CCSM)

            [var,zi,loni,lati] = get_CCSM_data(BGC_INI.bgc_frcini{BGC_INI.tracer_CCSM(trc)}, ...
                       BGC_INI.bgc_tracer{BGC_INI.tracer_CCSM(trc)},BRYtime.month(days),lon,lat,BRYtime,days);
            input_var_name{trc} = BGC_INI.bgc_tracer{BGC_INI.tracer_CCSM(trc)} ;

            dummy_mask = squeeze(var(end,:,:)) ; dummy_mask(dummy_mask~=0)=1 ;

               tic
               disp(' ')
               disp('Computing interpolation coefficients');
               [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati), ...
                          lon,lat,dummy_mask) ;
               A = get_hv_coef(zi, zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
               toc

           % fillmask
           var=fillmask(var, 1, dummy_mask, nnel);

           [nz,ny,nx] = size(var);

           %   Process scalar 3D variables
           var = reshape(A*reshape(var,nz*ny*nx,1),Nc,Mc,Lc);
           if (bnd ==1 )
              ncwrite(bryname,[input_var_name{trc} '_south'],permute(var(:, 1, :),[3 1 2]),[1 1 days])
           end
           if (bnd ==2 )
              ncwrite(bryname,[input_var_name{trc} '_east'],permute(var(:, :, end),[2 1 3]),[1 1 days])
           end
           if (bnd == 3 )
              ncwrite(bryname,[input_var_name{trc} '_north'],permute(var(:, end, :),[3 1 2]),[1 1 days])
           end
           if (bnd == 4 )
              ncwrite(bryname,[input_var_name{trc} '_west'],permute(var(:, :, 1),[2 1 3]),[1 1 days])
           end

      end    % End loop trc

  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (BGC_INI.data_glodapv2==1)

      clear vars vars_new input_var_name
      for trc=1:length(BGC_INI.tracer_glodapv2)

            [var,zi,loni,lati] = get_glodapv2_data(BGC_INI.bgc_frcini{BGC_INI.tracer_glodapv2(trc)}, ...
                       BGC_INI.bgc_tracer{BGC_INI.tracer_glodapv2(trc)},BRYtime.month(days),lon,lat,BRYtime,days);
            input_var_name{trc} = BGC_INI.bgc_tracer{BGC_INI.tracer_glodapv2(trc)} ;

            dummy_mask = squeeze(var(end,:,:)) ; dummy_mask(dummy_mask~=0)=1 ;

               tic
               disp(' ')
               disp('Computing interpolation coefficients');
               [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati), ...
                          lon,lat,dummy_mask) ;
               A = get_hv_coef(zi, zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
               toc

           % fillmask
           var=fillmask(var, 1, dummy_mask, nnel);

           [nz,ny,nx] = size(var);

           %   Process scalar 3D variables
           var = reshape(A*reshape(var,nz*ny*nx,1),Nc,Mc,Lc);
           if (bnd ==1 )
              ncwrite(bryname,[input_var_name{trc} '_south'],permute(var(:, 1, :),[3 1 2]),[1 1 days])
           end
           if (bnd ==2 )
              ncwrite(bryname,[input_var_name{trc} '_east'],permute(var(:, :, end),[2 1 3]),[1 1 days])
           end
           if (bnd == 3 )
              ncwrite(bryname,[input_var_name{trc} '_north'],permute(var(:, end, :),[3 1 2]),[1 1 days])
           end
           if (bnd == 4 )
              ncwrite(bryname,[input_var_name{trc} '_west'],permute(var(:, :, 1),[2 1 3]),[1 1 days])
           end

      end    % End loop trc

  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (BGC_INI.data_SYnn==1)

      clear vars vars_new input_var_name
      for trc=1:length(BGC_INI.tracer_SYnn)

            [var,zi,loni,lati] = get_SYnn_data(BGC_INI.bgc_frcini{BGC_INI.tracer_SYnn(trc)}, ...
                       BGC_INI.bgc_tracer{BGC_INI.tracer_SYnn(trc)},BRYtime.month(days),lon,lat,BRYtime,days);
            input_var_name{trc} = BGC_INI.bgc_tracer{BGC_INI.tracer_SYnn(trc)} ;

            dummy_mask = squeeze(var(end,:,:)) ; dummy_mask(dummy_mask~=0)=1 ;

               tic
               disp(' ')
               disp('Computing interpolation coefficients');
               [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati), ...
                          lon,lat,dummy_mask) ;
               A = get_hv_coef(zi, zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
               toc

%            % fillmask
%            var=fillmask(var, 1, dummy_mask, nnel);

           [nz,ny,nx] = size(var);

           %   Process scalar 3D variables
           var = reshape(A*reshape(var,nz*ny*nx,1),Nc,Mc,Lc);
           vars_new(:,:,:,trc) = var ;
      end    % End loop trc

            [C,trc1]=find(ismember(BGC_INI.bgc_tracer,'N2O_SIDEN')==1) ;
            if (C==1)
                [C_NO2 ,trc_NO2 ] = find(ismember(input_var_name,'N2O')==1) ;
                [C_NO2A,trc_NO2A] = find(ismember(input_var_name,'N2O_ATM')==1) ;
                if (C_NO2+C_NO2A<2)
                    error('ERROR computing N2O_SIDEN')
                else
                   N2O_SIDEN =  vars_new(:,:,:,trc_NO2)-vars_new(:,:,:,trc_NO2A) ;
                   N2O_SIDEN(N2O_SIDEN<0)=0;
                   vars_new(:,:,:,trc_NO2) = N2O_SIDEN + vars_new(:,:,:,trc_NO2A) ; % rescale N2O
%                    ncwrite(ininame,'N2O_SIDEN',permute(N2O_SIDEN,[3 2 1]),[icb jcb 1 1]) ;
                end
            end

      for trc=1:length(BGC_INI.tracer_SYnn)
           if (bnd ==1 )
              ncwrite(bryname,[input_var_name{trc} '_south'],permute(vars_new(:, 1, :,trc),[3 1 2]),[1 1 days])
           end
           if (bnd ==2 )
              ncwrite(bryname,[input_var_name{trc} '_east'],permute(vars_new(:, :, end,trc),[2 1 3]),[1 1 days])
           end
           if (bnd == 3 )
              ncwrite(bryname,[input_var_name{trc} '_north'],permute(vars_new(:, end, :,trc),[3 1 2]),[1 1 days])
           end
           if (bnd == 4 )
              ncwrite(bryname,[input_var_name{trc} '_west'],permute(vars_new(:, :, 1,trc),[2 1 3]),[1 1 days])
           end
      end    % End loop trc
      if (C==1)
      if (bnd ==1 )
         ncwrite(bryname,'N2O_SIDEN_south',permute(N2O_SIDEN(:, 1, :),[3 1 2]),[1 1 days])
      end
      if (bnd ==2 )
         ncwrite(bryname,'N2O_SIDEN_east',permute(N2O_SIDEN(:, :, end),[2 1 3]),[1 1 days])
      end
      if (bnd == 3 )
         ncwrite(bryname,'N2O_SIDEN_north',permute(N2O_SIDEN(:, end, :),[3 1 2]),[1 1 days])
      end
      if (bnd == 4 )
         ncwrite(bryname,'N2O_SIDEN_west',permute(N2O_SIDEN(:, :, 1),[2 1 3]),[1 1 days])
      end
      end

  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (BGC_INI.data_SeaWiFS==1)

      clear vars vars_new input_var_name
      for trc=1:length(BGC_INI.tracer_SeaWiFS)

            [var,loni,lati] = get_SeaWiFS_data(BGC_INI.bgc_frcini{BGC_INI.tracer_SeaWiFS(trc)}, ...
                       BGC_INI.bgc_tracer{BGC_INI.tracer_SeaWiFS(trc)},BRYtime.month(days),lon,lat,BRYtime,days);
            input_var_name{trc} = BGC_INI.bgc_tracer{BGC_INI.tracer_SeaWiFS(trc)} ;

            dummy_mask = squeeze(var(:,:)) ; dummy_mask(dummy_mask~=0)=1 ;

               tic
               disp(' ')
               disp('Computing interpolation coefficients');
               [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati), ...
                          lon,lat,dummy_mask) ;
               toc

%            % fillmask
%            var=fillmask(var, 1, dummy_mask, nnel);

%            [nz,ny,nx] = size(var);

           %   Process scalar 3D variables
           var(var<0.002)=0.002;
           var = var/0.75 ;
           var = sum(coef2d.*var(elem2d), 3) ;

           % construct vertical CHL profile
           disp('constructing vertical CHL profile from Morel and Berthon (1989)')
           vars_nnew = extr_chlo(var,zr) ;

           if (bnd ==1 )
              ncwrite(bryname,[input_var_name{trc} '_south'],permute(vars_nnew(:, 1, :),[3 1 2]),[1 1 days])
           end
           if (bnd ==2 )
              ncwrite(bryname,[input_var_name{trc} '_east'],permute(vars_nnew(:, :, end),[2 1 3]),[1 1 days])
           end
           if (bnd == 3 )
              ncwrite(bryname,[input_var_name{trc} '_north'],permute(vars_nnew(:, end, :),[3 1 2]),[1 1 days])
           end
           if (bnd == 4 )
              ncwrite(bryname,[input_var_name{trc} '_west'],permute(vars_nnew(:, :, 1),[2 1 3]),[1 1 days])
           end

      end    % End loop trc

  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (BGC_INI.data_Takahashi==1)

      clear vars vars_new input_var_name
      for trc=1:length(BGC_INI.tracer_Takahashi)

            [var,loni,lati] = get_Takahashi_data(BGC_INI.bgc_frcini{BGC_INI.tracer_Takahashi(trc)}, ...
                       BGC_INI.bgc_tracer{BGC_INI.tracer_Takahashi(trc)},BRYtime.month(days),lon,lat,BRYtime,days);
            var(:,end+1) = var(:,1);
            loni(:,end+1) = loni(:,1) + 360;
            lati(:,end+1) = lati(:,1);
            input_var_name{trc} = BGC_INI.bgc_tracer{BGC_INI.tracer_Takahashi(trc)} ;

            dummy_mask = squeeze(var(:,:)) ; dummy_mask(dummy_mask~=0)=1 ;

               tic
               disp(' ')
               disp('Computing interpolation coefficients');
               [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati), ...
                          lon,lat,dummy_mask) ;
               toc

%            % fillmask
%            var=fillmask(var, 1, dummy_mask, nnel);

%            [nz,ny,nx] = size(var);

           %   Process scalar 3D variables
           var = sum(coef2d.*var(elem2d), 3) ;

           if (bnd ==1 )
              ncwrite(bryname,[input_var_name{trc} '_south'],permute(var( 1, :),[2 1]),[1 days])
           end
           if (bnd ==2 )
              ncwrite(bryname,[input_var_name{trc} '_east'],permute(var( :, end),[1 2]),[1 days])
           end
           if (bnd == 3 )
              ncwrite(bryname,[input_var_name{trc} '_north'],permute(var(end, :),[2 1]),[1 days])
           end
           if (bnd == 4 )
              ncwrite(bryname,[input_var_name{trc} '_west'],permute(var( :, 1),[1 2]),[1 days])
           end

      end    % End loop trc

  end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [C,trc]=find(ismember(BGC_INI.bgc_tracer,'basindx')==1) ;
%  if (C==1)
%    [var] = get_sodazclm_data(BGC_INI.bgc_frcini{trc}, ...
%                                   BGC_INI.bgc_tracer{trc});
% %    ncwrite(ininame,BGC_INI.bgc_tracer{trc},permute(var,[2 1]),[1 1 1])
%    if (bnd ==1 )
%        ncwrite(bryname,[BGC_INI.bgc_tracer{trc} '_south'],permute(var( 1, :),[2 1]),[1 days])
%    end
%    if (bnd ==2 )
%        ncwrite(bryname,[BGC_INI.bgc_tracer{trc} '_east'],permute(var( :, end),[1 2]),[1 days])
%    end
%    if (bnd == 3 )
%        ncwrite(bryname,[BGC_INI.bgc_tracer{trc} '_north'],permute(var(end, :),[2 1]),[1 days])
%    end
%    if (bnd == 4 )
%        ncwrite(bryname,[BGC_INI.bgc_tracer{trc} '_west'],permute(var( :, 1),[1 2]),[1 days])
%    end
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end    % End loop bnd


  return























