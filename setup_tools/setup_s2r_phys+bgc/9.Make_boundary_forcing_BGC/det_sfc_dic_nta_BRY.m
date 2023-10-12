function [srf_NTA_vol,srf_dic_vol]=det_sfc_dic_nta_BRY(bryname,grdname,pars,bnd)

% dset_sfc_dic_nta.m.m: determine the surface boundary conditions for the carbon system
% based on pCO2 and other variables
% All operations are performed on the OA file, not the climatology file.
%
% The seawater package is used to calculate density.
%
% H. Frenzel, IGPP, UCLA
%
% First version (as DetCarbonBC.m): June 26, 2002
% Major revision (integrated into Roms_tools, eliminated all grid-specific
% settings etc.): December 21, 2007
%
% Revision: July 2009 Jean-Diego Santaren (JDS)
% The procedure can now deal with oa files used within make_carbon and with
% intermediary files used within make_biol_bec.
% Intermediate file has a coarser resolution than the 
% oa files and therefore soften the memory use
%

   if (bnd ==1 )
       suffix = '_south';
   end
   if (bnd ==2 )
      suffix = '_east';
   end
   if (bnd == 3 )
      suffix = '_north';
   end
   if (bnd == 4 )
      suffix = '_west';
   end 

%
% pCO2   D A T A
%
% get the pCO2 data
pco2 = ncread(bryname,['pCO2' suffix])'; 
pco2_siz = size(pco2);
L = pco2_siz(2);
M = pco2_siz(1);

%
% phosphate and silicate
%
srf_po4 = ncread(bryname,['PO4' suffix]) ; srf_po4 = squeeze(srf_po4(:,end,:))' ;
srf_sil = ncread(bryname,['SiO3' suffix]); srf_sil = squeeze(srf_sil(:,end,:))' ;

%
% T E M P   a n d   S A L T   d a t a
%
srf_temp = ncread(bryname,['temp' suffix]); srf_temp = squeeze(srf_temp(:,end,:))' ;
srf_salt = ncread(bryname,['salt' suffix]); srf_salt = squeeze(srf_salt(:,end,:))' ;


% Variables needed for the computations below
lon = ncread(grdname,'lon_rho')';
lon_mod360 = mod(lon,360);
lat = ncread(grdname,'lat_rho')';
[NX,NY] = size(lat);
NZ = pars.N ;

srf_density =  sw_dens0(srf_salt,srf_temp);

 if bnd==1 
   disp('South boundary')
   lon = lon(1,:) ;
   lat = lat(1,:) ;
   lon_mod360 = lon_mod360(1,:) ;
   for t=1:size(srf_density,1)
       lon2d(t,:) = lon ;
       lat2d(t,:) = lat ;
       lon2d_mod360(t,:) = lon_mod360 ;
   end
  end
  if bnd==2 
   disp('East boundary')
   lon = lon(:,end) ;
   lat = lat(:,end) ;
   lon_mod360 = lon_mod360(:,end) ;
   for t=1:size(srf_density,1)
       lon2d(t,:) = lon ;
       lat2d(t,:) = lat ;
       lon2d_mod360(t,:) = lon_mod360 ;
   end
  end
  if bnd==3 
   disp('North boundary')
   lon = lon(end,:) ;
   lat = lat(end,:) ;
   lon_mod360 = lon_mod360(end,:) ;
   for t=1:size(srf_density,1)
       lon2d(t,:) = lon ;
       lat2d(t,:) = lat ;
       lon2d_mod360(t,:) = lon_mod360 ;
   end
  end
  if bnd==4 
   disp('West boundary')
   lon = lon(:,1) ;
   lat = lat(:,1) ;
   lon_mod360 = lon_mod360(:,1) ;
   for t=1:size(srf_density,1)
       lon2d(t,:) = lon ;
       lat2d(t,:) = lat ;
       lon2d_mod360(t,:) = lon_mod360 ;
   end
  end


    % A L K A L I N I T Y   after Lee et al.: Global relationships of total alkalinity with salinity and temperature
    % in surface waters of the world's oceans, GRL, VOL. 33, L19605, doi:10.1029/2006GL027207, 2006)
    %
    % Table 1 distinguishes 5 regions,
    % each one has its own formula for NTA (umol kg-1)
    % read the basin index - Lee2006  doesn't have separate data for the
    % Mediterranean, so we'll have to use Atlantic data there as well
    
    disp('--- Surf Total Alkalinity using Lee et al. (2006) algorithm')
    
    basin_indx = ncread(bryname,['basindx' suffix])' ;
    
%     mask_atlantic = (basin_indx == 1 | basin_indx == 4);
%     mask_pacific = (basin_indx == 2);
%     mask_indian = (basin_indx == 3);
%     
%     % Zone 2 : Equatorial upwelling Pacific
%     %mm mask_zone2a = (lon_mod360 <= 285) .* (lon_mod360 >= 250) .* ...
%     %mm  (lat >= -20) .* (lat <= 20) .* mask_pacific;
%     %mm New we extend the equatorial Pac upwelling all the way to the west
%     %coast of America not stopping at 75W
%     mask_zone2a = (lon_mod360 >= 250) .* ...
%         (lat >= -20) .* (lat <= 20) .* mask_pacific;
%     mask_zone2b =  (lon_mod360 <= 250) .* (lon_mod360 >= 220) .* ...
%         (lat >= -10) .* (lat <= 10) .* mask_pacific;
%     mask_zone2 = mask_zone2a | mask_zone2b;
%     indx_zone2 = find(mask_zone2a | mask_zone2b);
%     num_2 = length(indx_zone2);
%     % Zone 1 : (Sub)tropics Atl, Ind, and Pac Oceancs excluding Zone 2
%     mask_zone1 = ( (lat >= -30) .* (lat <= 30) .* ...
%         (mask_atlantic + mask_indian + mask_pacific) .*...
%         ~mask_zone2a .* ~mask_zone2b ) ~= 0 ;
%     indx_zone1 = find(mask_zone1);
%     num_1 = length(indx_zone1);
%     % Zone 3 : North Atlantic
%     mask_zone3 = ((lat >= 30) .* mask_atlantic ) ~= 0;
%     indx_zone3 = find(mask_zone3);
%     num_3 = length(indx_zone3);
%     % Zone 4 : North Pacific
%     mask_zone4 = ( (lat >= 30) .* mask_pacific ) ~= 0;
%     indx_zone4 = find(mask_zone4);
%     num_4 = length(indx_zone4);
%     % Zone 5 : Southern Ocean
%     mask_zone5 = ( (lat <= -30) .* (mask_atlantic + mask_indian + mask_pacific) ) ~= 0;
%     indx_zone5 = find(mask_zone5);
%     num_5 = length(indx_zone5);
%     
%     sum_num = num_1 + num_2 + num_3 + num_4 + num_5;
%     sum_size = numel(basin_indx);
%     if (sum_num ~= sum_size)
%         disp(['WARNING in det_sfc_dic_nta: sum of indices in individual regions (', ...
%             num2str(sum_num), ') does not match total (', num2str(sum_size), ')']);
%         test = mask_zone1 + mask_zone2 + mask_zone3 + mask_zone4 + mask_zone5;
%         mask_zone4(test==0)=1 ;
%         indx_zone4 = find(mask_zone4);
%         num_4 = length(indx_zone4);
%     end
%     
%     % srf_NTA is in umol kg-1
%     srf_NTA=zeros(NX,NY);
%     TA     =zeros(NX,NY);
% 
%         srf_NTA(mask_zone1) = 2305 + 58.66 * (srf_salt(mask_zone1) - 35) + 2.32 * (srf_salt(mask_zone1) - 35).^2 - ...
%             1.41 * (srf_temp(mask_zone1) - 20) + 0.040 * (srf_temp(mask_zone1) - 20).^2;
%         %
%         srf_NTA(mask_zone2) = 2294 + 64.88 * (srf_salt(mask_zone2) - 35) + 0.39 * (srf_salt(mask_zone2) - 35).^2 - ...
%             4.52 * (srf_temp(mask_zone2) - 29) - 0.232 * (srf_temp(mask_zone2) - 29).^2;
%         %
%         srf_NTA(mask_zone3) = 2305 + 53.97 * (srf_salt(mask_zone3) - 35) + 2.74 * (srf_salt(mask_zone3) - 35).^2 - ...
%             1.16 * (srf_temp(mask_zone3) - 20) - 0.040 * (srf_temp(mask_zone3) - 20).^2;
%         %
%         %temp_lon=lon(indx_zone4);
%         
%         % modif IF&JDS
%         % Lee algo needs longitude to be between [0 360]
%         % and NOT between [-180 180]
%         %ii=find(temp_lon<0);
%         %if ~isempty(ii)
%         %    temp_lon(ii)=temp_lon(ii)+360;
%         %end
%         % end jds
%         %temp_lon = mod(temp_lon,360);
%         srf_NTA(mask_zone4) = 2305 + 53.23 * (srf_salt(mask_zone4) - 35) + 1.85 * (srf_salt(mask_zone4) - 35).^2 - ...
%             14.72 * (srf_temp(mask_zone4) - 20) - 0.158 * (srf_temp(mask_zone4) - 20).^2 + ...
%             0.062 * (srf_temp(mask_zone4) - 20) .* lon_mod360(mask_zone4);
%         srf_NTA(mask_zone5) = 2305 + 52.48 * (srf_salt(mask_zone5) - 35) + 2.85 * (srf_salt(mask_zone5) - 35).^2 - ...
%             0.49 * (srf_temp(mask_zone5) - 20) + 0.086 * (srf_temp(mask_zone5) - 20).^2;
%           
% %         % smooth the results
%              srf_NTA = lee_bdry_smooth(double(lon_mod360),double(lat), ...
%                              double(srf_NTA),double(basin_indx),double(2));
                         
        srf_NTA = lee_sfc_alk( double(srf_temp), ...
                               double(srf_salt), ...
                               double(lon2d), ...
                               double(lat2d), ...
                               double(basin_indx))  ;
                           
%         % smooth the results                   
%         srf_NTA = lee_bdry_smooth(double(lon2d_mod360),double(lat2d), ...
%                              double(srf_NTA),double(basin_indx),double(2));

%         
        % total (not normalized) alkalinity is needed for the calc_dic
        % routine
        TA = srf_NTA;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  IF: In Lee algorithm the surface alk is NOT normalized by salinity
        %     (contrary to Millero algo.)
        % we normalized it for subsequent computations
        
        srf_NTA = srf_NTA ./ srf_salt * 35.0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        srf_NTA_vol = srf_NTA ; % density is in kg/m3
%         srf_NTA_vol = 0.001 * srf_NTA .* srf_density; % density is in kg/m3
%         ncvarput(nc,'NTA_srf',srf_NTA_vol,[0,0,tt-1],[NX,NY,1]); % mmol m-3
        
%
% C R E A T E     D I C     C L I M A T O L O G Y
%
    disp('COMPUTING DIC  with calc_dic1')
    srf_dic = calc_dic1(double(TA), ...
                        double(pco2), ...
                        double(srf_temp), ...
                        double(srf_salt), ...
                        double(srf_po4), ...
                        double(srf_sil) );
    [NX,NY] = size(TA) ;              
    srf_dic = reshape(srf_dic,[NX,NY]);                

    % convert DIC from umol/kg to mmol/m3 (density is in kg/m3)
%     srf_dic_vol = 0.001 * srf_dic .* srf_density ;
    srf_dic_vol =  srf_dic ;
%     ncvarput(nc,'DIC_srf',srf_dic_vol,[0,0,tt-1],[NX,NY,1]); % mmol m-3

%


return
end
%
%
%
function srf_NTA = lee_bdry_smooth(lon,lat,srf_NTA,basin_indx,sband)
if nargin < 5
    sband = 2.5;
end
basin= mode(basin_indx(abs(lat)<=30));  % Basin is the most frequent basin_indx
% reset value to NaN  and use TriScatterInterp to refill these
% values
srf_NTA1 = srf_NTA;
srf_NTA(basin_indx ~= basin & abs(lat)<=30) = NaN; % Make sure we stick to our basin in the tropics
% This breaks if we set up a global model! check for that
srf_NTA(abs(abs(lat)-30)<=sband) = NaN; % 30N, 30S
srf_NTA(abs(abs(lat)-20)<=sband &...
    lon>=250 &... % extended to coast of SA lon<=285 &...
    basin_indx == 2) = NaN; % Pacific upwell bdr 1
srf_NTA(abs(abs(lat)-10)<=sband &...
    lon>=220 & lon<=250 & basin_indx == 2) = NaN; % Pacific upwell bdr 2
srf_NTA(abs(lon-250)<=sband & ...
    abs(lat)<=10 & basin_indx == 2 ) = NaN; % Pacific upwell bdr 2
% Upwelling extented to coast. Smoothing below no longer needed
%srf_NTA(abs(lon-285)<=sband & ...
%    abs(lat)<=20 & abs(lat)>=10& basin_indx == 2 ) = NaN; % Pacific upwell bdr 2
tofill = isnan(srf_NTA);
F = TriScatteredInterp(lon(~tofill),lat(~tofill),srf_NTA(~tofill),'natural');
srf_NTA(tofill) = F(lon(tofill),lat(tofill));
miss = isnan(srf_NTA);
if any(miss(:))
    F.Method = 'nearest';
    srf_NTA(miss) = F(lon(miss),lat(miss));
    % This needs fixing: 
    %%srf_NTA1 = filter2([2 3 2;3 4 3;2 3 2]/24.,srf_NTA1,'same');
    % The fix should be to interpolate linearly transition the different
    % parametrisations of Lee et al. AND GET RID of the smoothing in the
    % transition zones. Otherwise e.g. corners with minima will be filled with 
    % a high bias by any interpolation. This so happens in the south of the
    % USWC setup.  
    %%srf_NTA(miss) = srf_NTA1(miss);
end
return
end

