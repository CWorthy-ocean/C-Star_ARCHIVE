function [var,zi,loni,lati] = get_woa18_data(woa_file,variable,month,lon,lat,BRYtime,days) ;

% Processes and returns data from woa18 monthly climatology that is mask filled
% and extended from -180 to 360 degrees.

%%%% FIND the file to be read %%%%%%
   indstr = strfind(woa_file,'*') ;
   if month < 10
      month_str = ['0' num2str(month)];
   else
      month_str = num2str(month) ;
   end
   indir = dir([woa_file(1:indstr-1) month_str '*nc']);
   if (length(indir)>1)
      disp(['ERROR : cannot find the file to be read for ' variable]) ; error ;
   end
   infile = [indir.folder '/' indir.name];

  infile = woa_file ;
  disp(['get_woa18 : reading file : ' infile])

  indstr = strfind(woa_file,'annual');
  if isempty(indstr)==1
      indstr = strfind(woa_file,'seasonal');
  end
  infile_month = [woa_file(1:indstr-1) 'monthly_landfilled.nc'];
  disp(['get_woa18 : reading file : ' infile_month])

%%%% read lon/lat %%%%%%
  loni = ncread(infile,'lon')';
  disp('WARNING : special fix on longitude')
  loni = [loni(181:end) loni(1:180)+360];
  lati = ncread(infile,'lat')';
  i0 = find(loni<min(min(lon)),1,'last');
  i1 = find(loni>max(max(lon)),1,'first');
  j0 = find(lati<min(min(lat)),1,'last');
  j1 = find(lati>max(max(lat)),1,'first');
  [loni,lati] = meshgrid(loni,lati);
% figure(1);plot(lon,lat,'.k')
% figure(2);plot(loni,lati,'.k')
% error
  loni = loni(j0:j1,i0:i1);
  lati = lati(j0:j1,i0:i1);

  lev       = ncread(infile,'depth')';
%   lev_month = ncread(infile_month,'depth')';
  nz = length(lev);
% [month j0 j1 i0 i1]
  [ny,nx] = size(loni);
  zi = zeros(nz,ny,nx);
  for k = 1:nz
    zi(k,:,:) = -lev(k);
  end

  if strcmp(variable,'NO3')==1
      name='n_an'; factor = 1 ;
  elseif strcmp(variable,'PO4')==1
      name='p_an'; factor = 1 ;
  elseif strcmp(variable,'SiO3')==1
      name='i_an'; factor = 1 ;
  elseif strcmp(variable,'O2')==1
      name='o_an'; factor = 1000./ 22.4 ;
  else
      name=variable;
  end

  var_annual = ncread(infile,name) ;
  if size(var_annual,4)==1
     ind_an = [1 1 1 1 1 1 1 1 1 1 1 1];
  elseif size(var_annual,4)==4
     ind_an = [4 4 4 1 1 1 2 2 2 3 3 3];
  elseif size(var_annual,4)==12
     ind_an = [1 2 3 4 5 6 7 8 9 10 11 12];
  end
%%%  var_month  = ncread(infile_month,name) ;
  for t=1:12
      if strcmp(variable,'O2')
      var(:,:,:,t)    = var_annual(:,:,:,ind_an(t)) ;
%%%      var(:,:,1:57,t) = var_month(:,:,1:57,t)  ;
      else
      var(:,:,:,t)    = var_annual(:,:,:,ind_an(t)) ;
%%%      var(:,:,1:37,t) = var_month(:,:,1:37,t)  ;
      end
  end
  if length(size(var))==3
    var = permute(squeeze(var),[3 2 1]);
    var = permute([permute(var(:,:,181:end),[1,3,2]) permute(var(:,:,1:180),[1,3,2])],[1,3,2]);
    var = var(:,j0:j1,i0:i1)*factor;
  elseif length(size(var))==4
    if ( (month == 1) && BRYtime.day(days)<=(BRYtime.mdays(month)/2) )
       var = permute( BRYtime.interp_ratio(days) *squeeze(var(:,:,:,month)) + ...
                   (1-BRYtime.interp_ratio(days))*squeeze(var(:,:,:,12   )),[3 2 1]);
    elseif ( (month == 1) && BRYtime.day(days)>(BRYtime.mdays(month)/2) )
       var = permute( BRYtime.interp_ratio(days) *squeeze(var(:,:,:,month)) + ...
                   (1-BRYtime.interp_ratio(days))*squeeze(var(:,:,:,2    )),[3 2 1]);
    elseif ( (month == 12) && BRYtime.day(days)<=(BRYtime.mdays(month)/2) )
       var = permute( BRYtime.interp_ratio(days) *squeeze(var(:,:,:,month)) + ...
                   (1-BRYtime.interp_ratio(days))*squeeze(var(:,:,:,11   )),[3 2 1]);
    elseif ( (month == 12) && BRYtime.day(days)>(BRYtime.mdays(month)/2) )
       var = permute( BRYtime.interp_ratio(days) *squeeze(var(:,:,:,month)) + ...
                   (1-BRYtime.interp_ratio(days))*squeeze(var(:,:,:,1    )),[3 2 1]);
    elseif ( BRYtime.day(days)<=(BRYtime.mdays(month)/2) )
       var = permute( BRYtime.interp_ratio(days) *squeeze(var(:,:,:,month)) + ...
                   (1-BRYtime.interp_ratio(days))*squeeze(var(:,:,:,month-1)),[3 2 1]);
    elseif ( BRYtime.day(days)>(BRYtime.mdays(month)/2) )
       var = permute( BRYtime.interp_ratio(days) *squeeze(var(:,:,:,month)) + ...
                   (1-BRYtime.interp_ratio(days))*squeeze(var(:,:,:,month+1)),[3 2 1]);
    end
    var = permute([permute(var(:,:,181:end),[1,3,2]) permute(var(:,:,1:180),[1,3,2])],[1,3,2]);
    var = var(:,j0:j1,i0:i1)*factor;
  else
    disp('error reading WOA file') ; error ;
  end

%   disp('WARNING : fix for woa : bottom')
%   for i=1:i1-i0+1
%   for j=1:j1-j0+1
%       indnan = min(find(isnan(squeeze(var(:,j,i))))) ;
%       if (indnan~=1)
%       var(indnan:nz,j,i) = var(indnan-1,j,i);
%       end
%   end
%   end
%   var(isnan(var))=0;

disp('WARNING : return vertical in the ascending order')
for k=1:nz
    var_new (k,:,:) = var (nz-k+1,:,:) ;
    zi_new  (k,:,:) = zi  (nz-k+1,:,:) ;
end
var = var_new ;
zi  = zi_new  ;

end

% mypcolor(loni,lati,squeeze(u(end,:,:)));set(gca,'ydir','normal');colorbar
% mypcolor(loni,lati,ssh);set(gca,'ydir','normal');colorbar
% load /batavia/nmolem/OBSERV/COAST/coast_atl
% hold on;plot(clon,clat,'k');hold off
% error('testing get soda data')



