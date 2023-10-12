   tname = '/glade/scratch/bachman/ROMS_tools/DATASETS/tpxo9.v2a.nc';
   gname = '/glade/scratch/bachman/ROMS_tools/setup_r2r_phys+bgc/1.Make_grid/Iceland1_grd.nc';
   fname = 'Iceland1_tides.nc';

   nc = 10; % number of tidal constituents

 dn0 = datenum(1992,1,1);

 % Reference time of simulation
 year = 2000;
 month=    1;
 day  =    1;
 dn = datenum(year,month,day);

date_mjd=mjd(year,month,day);
[pf,pu,t0,phase_mkB]=egbert_correc(date_mjd,0,0,0);
pu = pu*pi/180;
phase_mkB = phase_mkB*pi/180;
aa = phase_mkB;

t = (dn-dn0)*3600*24;


lon = ncread(gname,'lon_rho');
lat = ncread(gname,'lat_rho');
  h = ncread(gname,'h');
mask= ncread(gname,'mask_rho');
ang = ncread(gname,'angle');
[nx,ny] = size(lon);

 % create forcing file
 om = ncread(tname,'omega',[1],[nc]);
 con= ncreadatt(tname,'/','Constituents');
 if ~exist(fname)
   create_frc_tides(gname,fname,nc,om,con)
   ncwriteatt(fname,'/','Reference_time',datestr(dn));
 end

%  find limits in the tpxo data file
 lon(lon<0) = lon(lon<0)+360;  % tpxo is between 0 and 360;

 lon0 = min(lon(:));
 lon1 = max(lon(:));
 lat0 = min(lat(:));
 lat1 = max(lat(:));

 tlonr = ncread(tname,'lon_r',[1 1],[inf 1]);
 tlatr = ncread(tname,'lat_r',[1 1],[1 inf]);

 % u,v, and r coordinates in the tpxo file are slightly shifted
 % the plus 1 is a bit of a hack, better would be to have different
 % i0,i1, etc for each.
 i0 = find(tlonr<lon0,1,'last');
 i1 = find(tlonr>lon1,1,'first')+1;
 j0 = find(tlatr<lat0,1,'last');
 j1 = find(tlatr>lat1,1,'first')+1;
 tnx = i1-i0+1;
 tny = j1-j0+1;

 tlonr = ncread(tname,'lon_r',[i0 j0],[tnx 1]);
 tlonu = ncread(tname,'lon_u',[i0 j0],[tnx 1]);
 tlonv = ncread(tname,'lon_v',[i0 j0],[tnx 1]);
 tlatr = ncread(tname,'lat_r',[i0 j0],[1 tny]);
 tlatu = ncread(tname,'lat_u',[i0 j0],[1 tny]);
 tlatv = ncread(tname,'lat_v',[i0 j0],[1 tny]);

 tlon2d = ncread(tname,'lon_r',[i0 j0],[tnx tny]);
 tlat2d = ncread(tname,'lat_r',[i0 j0],[tnx tny]);

 % Read elevations,sal, and transports from tpxo file
 thr = ncread(tname,'h_Re',[i0 j0 1],[tnx tny nc]);
 thi = ncread(tname,'h_Im',[i0 j0 1],[tnx tny nc]);
 thc = complex(thr,thi);
 tur = ncread(tname,'u_Re',[i0 j0 1],[tnx tny nc]);
 tui = ncread(tname,'u_Im',[i0 j0 1],[tnx tny nc]);
 tuc = complex(tur,tui);
 tvr = ncread(tname,'v_Re',[i0 j0 1],[tnx tny nc]);
 tvi = ncread(tname,'v_Im',[i0 j0 1],[tnx tny nc]);
 tvc = complex(tvr,tvi);
 tsr = ncread(tname,'sal_Re',[i0 j0 1],[tnx tny nc]);
 tsi = ncread(tname,'sal_Im',[i0 j0 1],[tnx tny nc]);
 Afact = 2.0;
 tsc = Afact*complex(tsr,tsi);   %% The Alan factor is 1.5-2.0
 clear 'thr' 'thi' 'tur' 'tui' 'tvr' 'tvi' 'tsr' 'tsi'

 % Get equilibrium tides and correct for SAL
 tpc = equi_tide(tlon2d,tlat2d,nc);
 tpc = tpc - tsc;

if 0
 figure(2)
 subplot(2,1,1)
 pcolor(tlon2d-360,tlat2d,abs(tpc(:,:,1)));shading flat;colorbar;colormap(jet)
 hold on
 plot(clon,clat,'k')
 plot(clon1,clat1,'k')
 hold off
 xlim([-180 -50])
 ylim([-50 60])
 subplot(2,1,2)
 pcolor(tlon2d-360,tlat2d,-180/pi*angle(tpc(:,:,1)));shading flat;colorbar;colormap(jet)
 hold on
 plot(clon,clat,'k')
 plot(clon1,clat1,'k')
 hold off
 xlim([-180 -50])
 ylim([-50 60])
 return
end

 % Multiply with nodal corrections and phase shifts to the reference time
 cI = complex(0,1);
 for ic = 1:nc
  thc(:,:,ic) = pf(ic)*thc(:,:,ic)*exp(cI*(om(ic)*t + pu(ic) + aa(ic)));
  tuc(:,:,ic) = pf(ic)*tuc(:,:,ic)*exp(cI*(om(ic)*t + pu(ic) + aa(ic)));
  tvc(:,:,ic) = pf(ic)*tvc(:,:,ic)*exp(cI*(om(ic)*t + pu(ic) + aa(ic)));
  tpc(:,:,ic) = pf(ic)*tpc(:,:,ic)*exp(cI*(om(ic)*t + pu(ic) + aa(ic)));
 end

 if 0
   % test time series for San Fransisco bar: 37.76°N  122.8°W   ->   237.2 E
   % https://tidesandcurrents.noaa.gov/noaatidepredictions.html?id=9414290
   % Match the dates to start at the reference time for comparison
   % Use meters and GMT time
   i = 946;
   j = 539;
   [tlonr(i) tlatr(j)]
   sshc = thc(i,j,:);
   tim = [0:0.01:2]*24*3600;
   ssh = 0*tim;
   for ic = 1:nc
    ssh = ssh + real( sshc(ic)*exp( cI*(om(ic)*tim)) );
   end
   plot(ssh+1)
   return
 end

 % Process variables and write to forcing file
 hc = zeros(nx,ny,nc);
 for ic = 1:nc
   hc(:,:,ic) = interp2(tlonr,tlatr,thc(:,:,ic).',lon,lat);  % non-conjugate transpose: .'
 end
 ncwrite(fname,'ssh_Re',real(hc));
 ncwrite(fname,'ssh_Im',imag(hc));

 % Process tidal barotropic fluxes
 uc = zeros(nx,ny,nc);
 vc = zeros(nx,ny,nc);
 for ic = 1:nc
   uc(:,:,ic) = interp2(tlonu,tlatu,tuc(:,:,ic).',lon,lat);
   vc(:,:,ic) = interp2(tlonv,tlatv,tvc(:,:,ic).',lon,lat);
 end
 clear 'tuc' 'tvc'

 % Rotate to grid orientation and convert to barotropic velocity
 cosa = cos(ang);
 sina = sin(ang);
 u = zeros(nx,ny,nc);
 v = zeros(nx,ny,nc);
 for ic = 1:nc
   u(:,:,ic) = (uc(:,:,ic).*cosa + vc(:,:,ic).*sina)./h;
   v(:,:,ic) = (vc(:,:,ic).*cosa - uc(:,:,ic).*sina)./h;
 end

 % Average to u and v points and write to file
 u = 0.5*(u(2:end,:,:)+u(1:end-1,:,:));
 v = 0.5*(v(:,2:end,:)+v(:,1:end-1,:));
 ncwrite(fname,'u_Re',real(u));
 ncwrite(fname,'u_Im',imag(u));
 ncwrite(fname,'v_Re',real(v));
 ncwrite(fname,'v_Im',imag(v));

 pc = zeros(nx,ny,nc);
 for ic = 1:nc
   pc(:,:,ic) = interp2(tlonr,tlatr,tpc(:,:,ic).',lon,lat);
 end
 ncwrite(fname,'pot_Re',real(pc));
 ncwrite(fname,'pot_Im',imag(pc));
 ncwriteatt(fname,'/','Alan factor of Tidal Potential',Afact);

return

