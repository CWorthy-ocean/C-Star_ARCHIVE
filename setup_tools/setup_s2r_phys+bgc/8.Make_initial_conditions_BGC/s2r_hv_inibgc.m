function s2r_hv_inibgc(BGC_INI,grdname,ininame,chdscd);
%--------------------------------------------------------------
%  Produce a ROMS initial file from BGC files data for January
%
%  Inspired by Roms_tools (IRD).
%  Thanks to Pierrick, Patrick (IRD), Francois (UCLA), Yusuke (UCLA)
%  Jeroen Molemaker (UCLA); nmolem@ucla.edu
%  adapted by pierre damien (ucla)
%--------------------------------------------------------------

BRYtime.day = 1 ;
BRYtime.interp_ratio = 1 ;
BRYtime.mdays = [31 28 31 30 31 30 31 31 30 31 30 31] ;

% Get S-coordinate params for child grid
  theta_b = chdscd.theta_b;
  theta_s = chdscd.theta_s;
  hc      = chdscd.hc;
  Nc      = chdscd.N;
  scoord_c  = chdscd.scoord;
  theta_s = ncread(ininame,'theta_s');
  theta_b = ncread(ininame,'theta_b');
  hc      = ncread(ininame,'hc');

  [npc mpc] = size(ncread(grdname,'h'));

  disp('WARNING : BGC initialisation set to january')
  month = 1;
%   time  = 15*3600*24;
  time  = 0;

% Set bry_time

  % Get child grid and chunk size
    ndomx = 1;
    ndomy = 1;

    [Lp,Mp] = size(ncread(grdname,'h'));

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
        h    = ncread(grdname,'h'       ,[icb jcb],[ice-icb+1 jce-jcb+1])';
        mask = ncread(grdname,'mask_rho',[icb jcb],[ice-icb+1 jce-jcb+1])';
        lon  = ncread(grdname,'lon_rho' ,[icb jcb],[ice-icb+1 jce-jcb+1])';
        lat  = ncread(grdname,'lat_rho' ,[icb jcb],[ice-icb+1 jce-jcb+1])';
        lon(lon<0) = lon(lon<0) + 360;

        [Mc,Lc] = size(mask);
        maskc3d = zeros(Nc,Mc,Lc);
        for k = 1:Nc
         maskc3d(k,:,:) = mask;
        end

        % Z-coordinate (3D) on child grid
        zr = zlevs4(h, h*0, theta_s, theta_b, hc, Nc, 'r', scoord_c);
        zw = zlevs4(h, h*0, theta_s, theta_b, hc, Nc, 'w', scoord_c);

        if (BGC_INI.data_woa18==1)

            clear vars vars_new input_var_name
            for trc=1:length(BGC_INI.tracer_woa18)

%             %Special traitment for O2 (different depth verctor)
%             if (strcmp(BGC_INI.bgc_tracer{BGC_INI.tracer_woa18(trc)},'O2')==1)

            [var,zi,loni,lati] = get_woa18_data(BGC_INI.bgc_frcini{BGC_INI.tracer_woa18(trc)}, ...
                                  BGC_INI.bgc_tracer{BGC_INI.tracer_woa18(trc)},month,lon,lat,BRYtime,1);
            input_var_name{trc} = BGC_INI.bgc_tracer{BGC_INI.tracer_woa18(trc)} ;

            tic
            disp('Computing interpolation coefficients');
            dummy_mask = squeeze(var(1,:,:))*0+1 ;
            [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
            A = get_hv_coef( zi, zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
            toc

%             % fillmask
%             var=fillmask(var, 1, dummy_mask, nnel);
            % interp
            [nz,ny,nx] = size(var);
            var = reshape(A*reshape(var,nz*ny*nx,1),Nc,Mc,Lc);

            disp(' Writing ini file')
            ncwrite(ininame,input_var_name{trc},permute(var ,[3 2 1]),[icb jcb 1 1])

            end

        end
        if (BGC_INI.data_CCSM==1)

            clear vars vars_new input_var_name
            for trc=1:length(BGC_INI.tracer_CCSM)
            [vars(:,:,:,trc),zi,loni,lati] = get_CCSM_data(BGC_INI.bgc_frcini{BGC_INI.tracer_CCSM(trc)}, ...
                                  BGC_INI.bgc_tracer{BGC_INI.tracer_CCSM(trc)},month,lon,lat,BRYtime,1);
            input_var_name{trc} =  BGC_INI.bgc_tracer{BGC_INI.tracer_CCSM(trc)} ;
            end

            tic
            disp('Computing interpolation coefficients');
            dummy_mask = squeeze(vars(1,:,:,1)) ; dummy_mask(dummy_mask~=0)=1 ;
            [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
            A = get_hv_coef( zi, zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
            toc
            % fillmask
            for trc=1:length(input_var_name)
              vars(:,:,:,trc)= fillmask(squeeze(vars(:,:,:,trc)), 1, dummy_mask, nnel);
            end
            % interp
            [nz,ny,nx,ntrc] = size(vars);
            for trc=1:length(input_var_name)
              vars_new(:,:,:,trc)= reshape(A*reshape(squeeze(vars(:,:,:,trc)),nz*ny*nx,1),Nc,Mc,Lc);
            end

            disp(' Writing ini file')

            for trc=1:length(input_var_name)
	        ncwrite(ininame,input_var_name{trc},permute(squeeze(vars_new(:,:,:,trc)),[3 2 1]),[icb jcb 1 1])
            end

        end
        if (BGC_INI.data_glodapv2==1)

            clear vars vars_new input_var_name
            for trc=1:length(BGC_INI.tracer_glodapv2)
            [vars(:,:,:,trc),zi,loni,lati] = get_glodapv2_data(BGC_INI.bgc_frcini{BGC_INI.tracer_glodapv2(trc)}, ...
                                  BGC_INI.bgc_tracer{BGC_INI.tracer_glodapv2(trc)},month,lon,lat,BRYtime,1);
            input_var_name{trc} =  BGC_INI.bgc_tracer{BGC_INI.tracer_glodapv2(trc)} ;
            end

            tic
            disp('Computing interpolation coefficients');
            dummy_mask = squeeze(vars(1,:,:,1)) ; dummy_mask(dummy_mask~=0)=1 ;
            [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
            A = get_hv_coef( zi, zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
            toc
            % fillmask
            for trc=1:length(input_var_name)
              vars(:,:,:,trc)= fillmask(squeeze(vars(:,:,:,trc)), 1, dummy_mask, nnel);
            end
            % interp
            [nz,ny,nx,ntrc] = size(vars);
            for trc=1:length(input_var_name)
              vars_new(:,:,:,trc)= reshape(A*reshape(squeeze(vars(:,:,:,trc)),nz*ny*nx,1),Nc,Mc,Lc);
            end

            disp(' Writing ini file')

            for trc=1:length(input_var_name)
	        ncwrite(ininame,input_var_name{trc},permute(squeeze(vars_new(:,:,:,trc)),[3 2 1]),[icb jcb 1 1])
            end

        end
        if (BGC_INI.data_SYnn==1)

            clear vars vars_new input_var_name
            for trc=1:length(BGC_INI.tracer_SYnn)
            [vars(:,:,:,trc),zi,loni,lati] = get_SYnn_data(BGC_INI.bgc_frcini{BGC_INI.tracer_SYnn(trc)}, ...
                                  BGC_INI.bgc_tracer{BGC_INI.tracer_SYnn(trc)},month,lon,lat,BRYtime,1);
            input_var_name{trc} =  BGC_INI.bgc_tracer{BGC_INI.tracer_SYnn(trc)} ;
            end

            tic
            disp('Computing interpolation coefficients');
            dummy_mask = squeeze(vars(1,:,:,1))*0 + 1 ;
            [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
            A = get_hv_coef( double(zi), zr, coef2d, elem2d, double(loni), double(lati), lon, lat);
            toc
            % interp
            [nz,ny,nx,ntrc] = size(vars);
            for trc=1:length(input_var_name)
              vars_new(:,:,:,trc)= reshape(A*reshape(squeeze(vars(:,:,:,trc)),nz*ny*nx,1),Nc,Mc,Lc);
            end

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
                   ncwrite(ininame,'N2O_SIDEN',permute(N2O_SIDEN,[3 2 1]),[icb jcb 1 1]) ;
                end
            end

            disp(' Writing ini file')

            for trc=1:length(input_var_name)
	        ncwrite(ininame,input_var_name{trc},permute(squeeze(vars_new(:,:,:,trc)),[3 2 1]),[icb jcb 1 1])
            end

        end
        if (BGC_INI.data_SeaWiFS==1)

            clear vars vars_new input_var_name vars_nnew
            for trc=1:length(BGC_INI.tracer_SeaWiFS)
            [vars(:,:,trc),loni,lati] = get_SeaWiFS_data(BGC_INI.bgc_frcini{BGC_INI.tracer_SeaWiFS(trc)}, ...
                                  BGC_INI.bgc_tracer{BGC_INI.tracer_SeaWiFS(trc)},month,lon,lat,BRYtime,1);
            input_var_name{trc} =  BGC_INI.bgc_tracer{BGC_INI.tracer_SeaWiFS(trc)} ;
            end

            tic
            disp('Computing interpolation coefficients');
            dummy_mask = squeeze(vars(:,:,1))*0 + 1 ;
            [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
            toc

            % interp
            for trc=1:length(input_var_name)
              var_test = squeeze(vars(:,:,trc)) ;
              var_test(var_test<0.002)=0.002;
              var_test = var_test/0.75 ;
              vars_new(:,:,trc) = sum(coef2d.*var_test(elem2d), 3) ;
            end

            % construct vertical CHL profile
            disp('constructing vertical CHL profile from Morel and Berthon (1989)')
            for trc=1:length(input_var_name)
              var_test = squeeze(vars_new(:,:,trc)) ;
              vars_nnew(:,:,:,trc) = extr_chlo(var_test,zr) ;
            end

            disp(' Writing ini file')

            for trc=1:length(input_var_name)
	        ncwrite(ininame,input_var_name{trc},permute(squeeze(vars_nnew(:,:,:,trc)),[3 2 1]),[icb jcb 1 1])
            end

        end
        if (BGC_INI.data_Takahashi==1)

            clear vars vars_new input_var_name vars_nnew
            for trc=1:length(BGC_INI.tracer_Takahashi)
            [var,loni,lati] = get_Takahashi_data(BGC_INI.bgc_frcini{BGC_INI.tracer_Takahashi(trc)}, ...
                                  BGC_INI.bgc_tracer{BGC_INI.tracer_Takahashi(trc)},month,lon,lat,BRYtime,1);
            input_var_name{trc} =  BGC_INI.bgc_tracer{BGC_INI.tracer_Takahashi(trc)} ;
            tic
            disp('Computing interpolation coefficients');
            dummy_mask = var*0 + 1 ;
            [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
            toc
            var = sum(coef2d.*var(elem2d), 3) ;
            ncwrite(ininame,input_var_name{trc},permute(var,[2 1]),[icb jcb 1])
            end

        end

      end
    end

   [C,trc]=find(ismember(BGC_INI.bgc_tracer,'basindx')==1) ;
   if (C==1)
   disp('We DONT NEED BASINDX SET TO 0')
%   [var,loni,lati] = get_sodazclm_data(BGC_INI.bgc_frcini{trc}, ...
%                                  BGC_INI.bgc_tracer{trc});
%   tic
%   disp('Computing interpolation coefficients');
%   dummy_mask = squeeze(var(:,:))*0 + 1 ;
   lon  = ncread(grdname,'lon_rho')';
   lat  = ncread(grdname,'lat_rho')';
   lon(lon<0) = lon(lon<0) + 360;
%   [elem2d,coef2d,nnel] = get_tri_coef(double(loni),double(lati),lon,lat,dummy_mask);
%   toc
%   var = sum(coef2d.*var(elem2d), 3) ;
%   ncwrite(ininame,BGC_INI.bgc_tracer{trc},permute(var,[2 1]),[1 1 1])
    ncwrite(ininame,BGC_INI.bgc_tracer{trc},permute(lon*0,[2 1]),[1 1 1])
   end

%   [C,trc]=find(ismember(BGC_INI.bgc_tracer,'DIC')==1) ;
%   if (C==1)
%       disp('Computing DIC and Alk from surface pCO2')
%       %[srf_NTA_vol,srf_dic_vol]=det_sfc_dic_nta(ininame,grdname,chdscd) ;
%       dic_alk_profile(ininame,grdname,chdscd,srf_NTA_vol,srf_dic_vol,-750) ;
%   end

   vec = strfind(BGC_INI.bgc_frctype,'s') ;
   for trc=1:length(vec)
       if vec{trc}==1
       scaled_tracer = BGC_INI.bgc_tracer{trc} ;
       ref_ind       = str2num(BGC_INI.bgc_frctype{trc}(3:end)) ;
       tracer_scaled = BGC_INI.bgc_tracer{ref_ind} ;
       factor        = str2num(BGC_INI.bgc_frcini{trc}) ;
       disp(['SCALING variable ' scaled_tracer ' from ' tracer_scaled ' and factor ' num2str(factor)])
       var = ncread(ininame,tracer_scaled) ;
       ncwrite(ininame,scaled_tracer,var*factor,[1 1 1 1]);
       end
   end


return
