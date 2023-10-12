function dic_alk_profile(ininame,grdname,pars,srf_NTA_vol,srf_dic_vol,depth_seas)
% DIC_ALK_PROFILE: GLODAP for annual, deeper layer
%  above depths_seas: Linear interpolation in density space from GLODAP
%  to Lee PCO2 based estimates of surface Alk/DIC
% USAGE
%    dic_alk_profile(grdfile,oafile,oabiofile,depth_seas[,dic_offset])
%   INPUT:
%       grdfile: grid file for mask_rho masking
%       oafile: temp/salt OA file  for density
%       depth_seas: depth [m] up to which to extent the monthly data (seasonality layer)
%       dic_offset: optional surface DIC offset (default: 0)
%
%
%  Replaces  det_dic_profile/ det_alk_profile routines by H. Frenzel, IGPP, UCLA
%
% M. Munnich (2011)
%
%

mask = ncread(grdname,'mask_rho')';
h    = ncread(grdname,'h')';
[NY,NX]=size(mask);
NZ = pars.N ;

zr = zlevs4(h, h*0, pars.theta_s, pars.theta_b, pars.hc, pars.N, 'r', pars.scoord);

k_ann3D = zr*0;
k_ann3D(zr < depth_seas)=1;

%
% Write lower, annual-only levels
% a) DIC
DIC = permute(ncread(ininame,'DIC_glodap'),[3 2 1]); % annual data, (no seasonality)

% b) Same for Alk ("TA")
TA = permute(ncread(ininame,'Alk_glodap'),[3 2 1]);

clear data
salt0 = 35.0;
salt0i = 1./salt0;


% ncwrite(ininame,'DIC',permute(DIC*1.023 ,[3 2 1]),[1 1 1 1]);

    temp = permute(ncread(ininame,'temp'),[3 2 1]);
    salt = permute(ncread(ininame,'salt'),[3 2 1]);

    rho = sw_dens(salt,temp,-zr); %  z_grd needs temp/salt dimensions

    dataDIC = DIC ;
    dataAlk = TA ;

    disp('Computing DIC and Alk profiles')

    for i=1:NX
    for j=1:NY

    k_ann=find(k_ann3D(:,j,i)==0,1,'first')   ;

    if k_ann~=1
    clear temp

    drho = bsxfun(@minus,rho(k_ann:NZ,j,i),rho(k_ann-1,j,i));
    % alpha: density estimated fraction of climatological
    %        annual cycle strength at depth
    alpha = bsxfun(@rdivide,drho,drho(end));

    % Linear interpolation of DIC along monthly density variations
%     data_srf = ncvarget(nc,'DIC_srf',[0,0,tt-1],[NX,NY,1]);
%     data_srf = data_srf + dic_offset; % Giuliana sensitivity tapering
%     deltadata = data_srf - DIC(:,:,k_ann);
    deltadata = srf_dic_vol(j,i)- DIC(k_ann-1,j,i);
    dataDIC1       = bsxfun(@plus,DIC(k_ann-1,j,i), bsxfun(@times,alpha,deltadata));
    dataDIC(:,j,i) = (cat(1,DIC(1:k_ann-1,j,i),dataDIC1))*0.001.*rho(:,j,i);  % AF to volume concentrations [mmol m-3]

  % Same for NTA (not TA to take alk salinity dependence into account)
%     data_srf = ncvarget(nc,'NTA_srf',[0,0,tt-1],[NX,NY,1]);
    NTA_below = TA(k_ann-1,j,i)*salt0./salt(k_ann-1,j,i);
    deltadata = srf_NTA_vol(j,i)-NTA_below;
    dataAlk1       = bsxfun(@plus,NTA_below, bsxfun(@times,alpha,deltadata));
%     dataAlk1       = dataAlk1*salt0i.*salt(k_ann:NZ,j,i);           % back to TA
%     dataAlk(:,j,i) = cat(1,TA(1:k_ann-1,j,i),dataAlk1)*0.001.*rho(:,j,i);     % AF to volume concentrations [mmol m-3]
    dataAlk(:,j,i) = cat(1,TA(1:k_ann-1,j,i)*salt0./salt(k_ann-1,j,i),dataAlk1)*0.001.*rho(:,j,i).*salt0i.*salt(:,j,i);
    end

    end
    end

    for k = 1:NZ
      dataAlk(k,:,:) = inpaintn(squeeze(dataAlk(k,:,:)));
      dataDIC(k,:,:) = inpaintn(squeeze(dataDIC(k,:,:)));
    end

    ncwrite(ininame,'DIC',permute(dataDIC ,[3 2 1]),[1 1 1 1]);
    ncwrite(ininame,'Alk',permute(dataAlk ,[3 2 1]),[1 1 1 1]);


end



