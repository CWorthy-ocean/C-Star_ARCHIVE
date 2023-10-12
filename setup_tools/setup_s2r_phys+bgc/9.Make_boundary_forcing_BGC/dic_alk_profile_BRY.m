function dic_alk_profile_BRY(bryname,grdname,pars,srf_NTA_vol,srf_dic_vol,depth_seas,bnd)
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

zr1 = zlevs4(h, h*0, pars.theta_s, pars.theta_b, pars.hc, pars.N, 'r', pars.scoord);

   if (bnd ==1 )
       suffix = '_south';
       DIC = permute(ncread(bryname,['DIC_glodap' suffix]),[3 2 1]);
       TA = permute(ncread(bryname,['Alk_glodap' suffix]),[3 2 1]);
       for t=1:size(DIC,1)
           zr(t,:,:) = squeeze(zr1(:,1,:)) ;
       end
       DIC = permute(DIC,[2 1 3]) ;
       TA  = permute(TA ,[2 1 3]) ;
       zr  = permute(zr ,[2 1 3]) ;
   end
   if (bnd ==2 )
      suffix = '_east';
      DIC = permute(ncread(bryname,['DIC_glodap' suffix]),[3 2 1]);
      TA = permute(ncread(bryname,['Alk_glodap' suffix]),[3 2 1]);
      for t=1:size(DIC,1)
           zr(t,:,:) = squeeze(zr1(:,:,end)) ;
      end
      DIC = permute(DIC,[2 1 3]) ;
      TA  = permute(TA ,[2 1 3]) ;
      zr  = permute(zr ,[2 1 3]) ; 
   end
   if (bnd == 3 )
      suffix = '_north';
      DIC = permute(ncread(bryname,['DIC_glodap' suffix]),[3 2 1]);
      TA = permute(ncread(bryname,['Alk_glodap' suffix]),[3 2 1]);
      for t=1:size(DIC,1)
           zr(t,:,:) = squeeze(zr1(:,end,:)) ;
      end
      DIC = permute(DIC,[2 1 3]) ;
      TA  = permute(TA ,[2 1 3]) ;
      zr  = permute(zr ,[2 1 3]) ; 
   end
   if (bnd == 4 )
      suffix = '_west';
      DIC = permute(ncread(bryname,['DIC_glodap' suffix]),[3 2 1]);
      TA = permute(ncread(bryname,['Alk_glodap' suffix]),[3 2 1]);
      for t=1:size(DIC,1)
           zr(t,:,:) = squeeze(zr1(:,:,1)) ;
      end
      DIC = permute(DIC,[2 1 3]) ;
      TA  = permute(TA ,[2 1 3]) ;
      zr  = permute(zr ,[2 1 3]) ; 
   end 


k_ann3D = zr*0; 
k_ann3D(zr < depth_seas)=1; 

%
% Write lower, annual-only levels
% a) DIC
% DIC = permute(ncread(bryname,'DIC_glodap'),[3 2 1]); % annual data, (no seasonality)

% b) Same for Alk ("TA")
% TA = permute(ncread(bryname,'Alk_glodap'),[3 2 1]);

clear data 
salt0 = 35.0;
salt0i = 1./salt0;


% ncwrite(ininame,'DIC',permute(DIC*1.023 ,[3 2 1]),[1 1 1 1]);

    temp = permute(ncread(bryname,['temp' suffix]),[2 3 1]); 
    salt = permute(ncread(bryname,['salt' suffix]),[2 3 1]);
    
    rho = sw_dens(salt,temp,-zr); %  z_grd needs temp/salt dimensions
    
    dataDIC = DIC ;
    dataAlk = TA ;
    
    disp('Computing DIC and Alk profiles')
    
    for i=1:size(rho,3)
    for j=1:size(rho,2)
        
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
    
    if (bnd ==1 )
    ncwrite(bryname,'DIC_south',permute(dataDIC ,[3 1 2]),[1 1 1]);
    ncwrite(bryname,'Alk_south',permute(dataAlk ,[3 1 2]),[1 1 1]);
    end
    if (bnd ==2 )
    ncwrite(bryname,'DIC_east',permute(dataDIC ,[3 1 2]),[1 1 1]);
    ncwrite(bryname,'Alk_east',permute(dataAlk ,[3 1 2]),[1 1 1]);
    end
    if (bnd ==3 )
    ncwrite(bryname,'DIC_north',permute(dataDIC ,[3 1 2]),[1 1 1]);
    ncwrite(bryname,'Alk_north',permute(dataAlk ,[3 1 2]),[1 1 1]);
    end
    if (bnd ==4 )
    ncwrite(bryname,'DIC_west',permute(dataDIC ,[3 1 2]),[1 1 1]);
    ncwrite(bryname,'Alk_west',permute(dataAlk ,[3 1 2]),[1 1 1]);
    end

end



