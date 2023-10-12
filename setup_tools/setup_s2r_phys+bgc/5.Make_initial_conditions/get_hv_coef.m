function A = get_hv_coef(zp,zc,coef2d,elem2d,lonp,latp,lonc,latc)

% Function returns...
%
% Inputs:
%          interpolation coefficients from the 2d (horizontal) interpolation.
%          pointer array to the elements from which the 2d interpolation is computed
%          3d parent z-values (zp)
%          3d child z-values (zc)
%          Both zp and zc most be in ascending order !!!!! (for the k index)
%
% Outputs:  3d interpolation matrix A
%           use like: Fc = reshape(A*reshape(Fp,Np*Mp*LP,1),Nc,Mc,Lc)
%
%         (c) 2007 Jeroen Molemaker

%  Construct a sparse matrix A such that Fc = A*Fp; (with some reshapes)
%  Let A be the product of a horizontal interpolation Ah and a vertical
%  interpolation Av.

  gnomonic_option = 0 ;

  [Np Mp Lp] = size(zp);
  [Nc Mc Lc] = size(zc);
  ndimp = Np*Mp*Lp;
  ndimc = Nc*Mc*Lc;
  ndimt = (Np+1)*Mc*Lc;

%  Ah interpolates from a (np,mp,lp) grid to a (np+1,mc,lc) grid.
   ca = zeros(Np+1,Mc,Lc,3);
   ic = zeros(Np+1,Mc,Lc,3);
   jc = zeros(Np+1,Mc,Lc,3);

   disp('allocated for Ah arrays');

   sub2d = zeros(Mc,Lc,3);
   for i = 1:Mc
    for j = 1:Lc
      sub2d(i,j,:) = sub2ind([Mc Lc],i,j);
    end
   end

   for k = 2:Np+1
     ca(k,:,:,:) = coef2d;
     ic(k,:,:,:) = k + (Np+1)*(sub2d -1);
     jc(k,:,:,:) = (k-1) + Np*(elem2d-1);
   end

   width1 = max(max(lonp))-min(min(lonp));
   width2 = max(max(latp))-min(min(latp));
   width3 = max(max(lonc))-min(min(lonc));
   width4 = max(max(latc))-min(min(latc));
   width = max([width1 width2 width3 width4]);
   if ( width<100 && gnomonic_option )
    lon0 = mean(mean(mean(lonc)));
    lat0 = mean(mean(mean(latc)));
    [xp,yp] = gnomonic(lonp,latp,lon0,lat0);
    [xc,yc] = gnomonic(lonc,latc,lon0,lat0);
   else
    disp('no gno')
    xp = lonp; yp = latp;
    xc = lonc; yc = latc;
   end

   %% This is unfortunately a somewhat ugly fix to deal
   %% with gross differences in topography between the 2 grids.
   %% It adds an extra vertical point taken from the intermediate (horizontal)
   %% interpolation. So Ah works from (Lp,Mp,Np) to (Lp+1,Mc,Nc)
   %% We add an extra point below when the depth mismatch exceeds
   %% a criterium

%  zp(:,1,1) 
%  zc(:,1,1) 
   zp_low = squeeze(zp(1,:,:));
   zc_low = squeeze(zc(1,:,:));
   zt_low = sum(coef2d.*zp_low(elem2d),3);

   mismatch = 20000; %% When the child grid's lowest point at (i,j)
                    %% is more than 'mismatch' below the
                    %% interpolated parent point, we put an
                    %% extra point below it by means
                    %% of nearest neighbor horizontal extrapolation.

   hmin = max(max(zp_low));
   hmax = min(min(zp_low));
   nlev = round(abs(hmax/mismatch));
   tri_lev = cell(nlev,1);
   for lev = 1:nlev
     depth =  lev*hmax/(nlev+1);
     shallow = find(zp_low > depth);
     x_tmp = xp;
     y_tmp = yp;
     x_tmp(shallow)  = x_tmp(shallow) + 2;
     y_tmp(shallow)  = y_tmp(shallow) + 2;
     X_tmp(:,:,lev)  = [reshape(x_tmp,Mp*Lp,1) reshape(y_tmp,Mp*Lp,1) ];
     tri_tmp = delaunayn(X_tmp(:,:,lev));
     tri_lev{lev} = tri_tmp;
   end
   disp('--- fixing parent/child topo mismatch (may take a while...)')
   done = 10;
   for i = 1:Mc
    for j = 1:Lc
      if zt_low(i,j) > zc_low(i,j)+ mismatch  %% We need an extra point below zt_low
        lev = round(zc_low(i,j)*(nlev+1)/hmax);
        %display('--- sideways interpolation due to parent/child topo mismatch')
        lev = max(1,lev);
        lev = min(nlev,lev);
        nnel  = dsearchn(X_tmp(:,:,lev),tri_lev{lev},[xc(i,j) yc(i,j)]);

        ca(1,i,j,1) = 1.0;
        ic(1,i,j,:) = 1 + (Np+1)*(sub2d(i,j) -1);
        jc(1,i,j,:) = 1 + Np*(nnel-1);
      else                  %% Trivial extra point
        ca(1,i,j,:) = ca(2,i,j,:);
        ic(1,i,j,:) = 1 + (Np+1)*(sub2d(i,j) -1);
        jc(1,i,j,:) = jc(2,i,j,:);
      end
    end
    percent_done = round(100*i/Mc);
    if percent_done==done
      disp(['------ ' num2str(percent_done) '% done'])
      done = done + 10;
    end
   end
   %% End of adding a lower layer

   ic = reshape(ic,ndimt*3,1);
   jc = reshape(jc,ndimt*3,1);
   ca = reshape(ca,ndimt*3,1);
   Ah = sparse(ic,jc,ca,ndimt,ndimp);

   disp('done with Ah');

   zc_tmp = reshape(Ah*reshape(zp,ndimp,1),Np+1,Mc,Lc);
   zc_tmp(1,:,:) = zc_tmp(1,:,:) - 0.1;  %% Avoid double z point in every column.

   clear ic
   clear jc
   clear ca

   Av = get_v_coef(zc_tmp,zc);
   A  = Av*Ah;
   disp('--- all done!!')
   return
