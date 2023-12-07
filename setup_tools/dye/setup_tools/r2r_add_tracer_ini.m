
ininame = '/glade/scratch/bachman/UCLA-ROMS/Work/Iceland4_dye/INPUT/ens0/Iceland4_rst.20120608120000.nc';
grdname = '/glade/scratch/bachman/ROMS_tools/Iceland4_dye/EASY/Iceland4_grd.nc';

parscd.N       = 100 ;
parscd.theta_s = 5.0;
parscd.theta_b = 2.0;
parscd.hc      = 300 ;
parscd.scoord = 'new2012';


ncid = netcdf.open(ininame,'NC_WRITE');    % WRITE

s_rho_id = netcdf.inqDimID(ncid,'s_rho')
xi_rho_id = netcdf.inqDimID(ncid,'xi_rho')
eta_rho_id = netcdf.inqDimID(ncid,'eta_rho')
time_id = netcdf.inqDimID(ncid,'time')
%time_id = netcdf.inqDimID(ncid,'time')

[s_rho, s_rho_len] = netcdf.inqDim(ncid,s_rho_id)
[xi_rho, xi_rho_len] = netcdf.inqDim(ncid,xi_rho_id)
[eta_rho, eta_rho_len] = netcdf.inqDim(ncid,eta_rho_id)
[time, time_len] = netcdf.inqDim(ncid,time_id)

nccreate(ininame,'dye','Dimensions',{'xi_rho',xi_rho_len,'eta_rho', eta_rho_len, 's_rho',s_rho_len,'time',time_len},'datatype','double');
ncwriteatt(ininame,'dye','long_name','dye concentration');
ncwriteatt(ininame,'dye','units','nondimensional');

dye = zeros(xi_rho_len, eta_rho_len, s_rho_len, time_len);

%%%%%%%%% Create dye patch

h    = ncread(grdname,'h');

dilation_factor = 1.1;
vert_decay_scale = 5;

radius = 80;
amplitude = 0;
center = [1200, 580];

for i = (center(1) - radius):(center(1) + radius)
  for j = (center(2) - radius):(center(2) + radius)

    ind = 1;
    dye_depths(ind) = -1;
    while dye_depths(end) > -h(i,j)
      ind = ind + 1;
      dye_depths(ind) = dye_depths(ind-1)*dilation_factor;
    end
    dye_depths = flipud(dye_depths);
    %depth_function = exp(dye_depths / vert_decay_scale);
    depth_function = zeros(size(dye_depths));
    depth_function(dye_depths > -5) = 1;

    zr = zlevs4(h(i,j), h(i,j)*0, parscd.theta_s, parscd.theta_b, parscd.hc, parscd.N, 'r', parscd.scoord);

    for k = 1:length(dye_depths)
      dye_tmp(k) = amplitude * depth_function(k) * exp(-( (i - center(1))^2 + (j-center(2))^2)/radius);
    end

    for t = 1:time_len
      dye(i,j,:,t) =  interp1(dye_depths, dye_tmp, zr, 'linear', 0);
    end


    clear dye_depths
    clear dye_tmp

  end
end

ncwrite(ininame, 'dye', dye);

netcdf.close(ncid)
