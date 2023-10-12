 u_south = ncread('pachug_edata.nc','nwpac_south_u');
 v_south = ncread('pachug_edata.nc','nwpac_south_v');
 r_south = ncread('pachug_edata.nc','nwpac_south_r');

 nx = 3680;
 ny = 1920;
 pname = '/avatar/nmolem/PACHUG/pachug_grd.nc';
 cname = '/avatar/nmolem/NWPAC/nwpac_grd.nc';
 smsk = ncread(cname,'mask_rho');
 smsk = smsk(:,1);
 h = ncread(pname,'h');
 msk = ncread(pname,'mask_rho');
 vmsk = msk(:,2:end).*msk(:,1:end-1);
 ip = [0:nx+1];
 jp = [0:ny+1];
 [ip,jp] = meshgrid(ip,jp);ip=ip';jp=jp';

 ipv = [-0.5:nx+0.5];
 jpv = [0:ny];
 [ipv,jpv] = meshgrid(ipv,jpv);ipv=ipv';jpv=jpv';
 ipvm = ipv;
 jpvm = jpv;
 ipvm(vmsk>0) = [];
 jpvm(vmsk>0) = [];

 % these are in 'absolute' index space, the boundary starts at zero
 iv = v_south(:,1);
 jv = v_south(:,2);
 ivm = iv;
 jvm = jv;
 ivm(smsk>0) = [];
 jvm(smsk>0) = [];


   ivs = iv + 1.5;
   jvs = jv + 1.0;
   sms = ivs*0+ 99;
   for i = 1:length(ivs)
    ill = floor(ivs(i));
    jll = floor(jvs(i));
    sms(i) = sum(sum(vmsk(ill:ill+1,jll:jll+1)));
   end

  x = 1:length(sms);
  ivm = iv;
  jvm = jv;
  xm = x;
  ivm(sms==0) = [];
  jvm(sms==0) = [];
  xm(sms==0) = [];

   indx = find((sms-1).*smsk<0);
   iv(indx) = interp1(xm,ivm,x(indx),'nearest','extrap');
   jv(indx) = interp1(xm,jvm,x(indx),'nearest','extrap');
   for i = 1:length(indx)

     ivs(indx(i)) = interp1
   end


