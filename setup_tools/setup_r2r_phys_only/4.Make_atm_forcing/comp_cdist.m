function cdist=comp_cdist(grdname,disname,coarse_frc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Compute distance to nearest masked point
%  For use in fill_frc_runnoff and wind_oro
%
%  2010, Jeroen Molemaker (UCLA)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% USER-DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
%

 disp(' ')
 disp(' Read in the grid')
 if coarse_frc
  lon  = ncread(grdname,'lon_coarse');
  lat  = ncread(grdname,'lat_coarse');
  mask = ncread(grdname,'mask_coarse');
 else
  lon  = ncread(grdname,'lon_rho');
  lat  = ncread(grdname,'lat_rho');
  mask = ncread(grdname,'mask_rho');
 end
 [nx,ny] = size(lon);
 d2r = pi/180;


  cdist = 0*lon + 1e10;

  nch(1) = 1;
  nch(2) = 3;
  nch(3) = 4;
  for it = 1:1
    ncx = nch(it);
    ncy = nch(it);
    for ic = 1:ncx
      for jc = 1:ncy
        [ic jc]
        i0 = 1 + (ic-1)*ceil(nx/ncx);
        i1 = i0+ceil(nx/ncx) + 20;
        j0 = 1 + (jc-1)*ceil(ny/ncy);
        j1 = j0+ceil(ny/ncy) + 20;
        i1 = min(i1,nx); j1 = min(j1,ny);
        lons = lon(i0:i1,j0:j1);
        lats = lat(i0:i1,j0:j1);
        masks= mask(j0:i1,j0:j1);
	lab = 0*masks + 1; lab(2:end-1,2:end-1) = masks(1:end-2,2:end-1)+masks(3:end,2:end-1)+masks(2:end-1,1:end-2)+masks(2:end-1,3:end);
        mlon = lons;mlon(masks>0|lab<1) = [];
        mlat = lats;mlat(masks>0|lab<1) = [];
	disp('here')
        swt = sum(~isnan(mlon))>0

        if sum(~isnan(mlon))>0
         for j = j0:j1
		 [j j1]
          for i = i0:i1
            if mask(i,j) < 1 
              cdist(i,j) = 0;
            else
              dist =  gc_dist(lon(i,j)*d2r,lat(i,j)*d2r,mlon*d2r,mlat*d2r);
              mdist = min(min(dist));
              cdist(i,j) = min(mdist,cdist(i,j));
            end
          end
         end
        else
         cdist(i0:i1,j0:j1) = min(3e5,cdist(i0:i1,j0:j1));
   
        
        
        end
      end
    end
  end

  cdist(cdist>3e5) = 3e5;

  imagesc(cdist');axis xy;colorbar

  save(disname,'cdist');
  disp('computed distance file');
  return
