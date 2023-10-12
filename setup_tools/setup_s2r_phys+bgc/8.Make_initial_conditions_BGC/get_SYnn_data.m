function [var,zi,loni,lati] = get_SYnn_data(SYnn_file,variable,month,lon,lat,BRYtime,days);

% Processes and returns data from woa18 monthly climatology that is mask filled
% and extended from -180 to 360 degrees.

%%%% FIND the file to be read %%%%%%
  infile = SYnn_file ;
  disp(['get_SYnn : reading file : ' infile])  
  
%%%% read lon/lat %%%%%%
  if ( contains(infile,'n2ofromnn')==1 )     
  loni = ncread(infile,'lon')'; 
  disp('WARNING : special fix on longitude')
  loni = [loni(181:end) loni(1:180)+360];
  lati = ncread(infile,'lat')';
  elseif ( contains(infile,'fromWOA13')==1 ) 
  loni = ncread(infile,'LON')'; 
  disp('WARNING : special fix on longitude')
  loni = [loni(181:end) loni(1:180)+360];
  lati = ncread(infile,'LAT')';
  end 
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
  
  lev = ncread(infile,'DEPTH')';
  nz = length(lev);
% [month j0 j1 i0 i1]
  [ny,nx] = size(loni);
  zi = zeros(nz,ny,nx);
  for k = 1:nz
    zi(k,:,:) = -lev(k);
  end
  
  if strcmp(variable,'N2O_ATM')==1
      name='N2O_SATURATION_MMOLPM3';
  else
      name=variable;
  end
  
  var = permute(squeeze(ncread(infile,name)),[4 3 2 1]);
  vec = size(var) ;
%   var = squeeze(var(month,:,:,:)) ;
  if (vec(1) == 12)
  if ( (month == 1) && BRYtime.day(days)<=(BRYtime.mdays(month)/2) )
      var =  BRYtime.interp_ratio(days) *squeeze(var(month,:,:,:)) + ...
          (1-BRYtime.interp_ratio(days))*squeeze(var(   12,:,:,:)) ;
  elseif ( (month == 1) && BRYtime.day(days)>(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month,:,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(    2,:,:,:)) ;           
  elseif ( (month == 12) && BRYtime.day(days)<=(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month,:,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(   11,:,:,:)) ; 
  elseif ( (month == 12) && BRYtime.day(days)>(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month,:,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(    1,:,:,:)) ;            
  elseif ( BRYtime.day(days)<=(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month  ,:,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(month-1,:,:,:)) ;  
  elseif ( BRYtime.day(days)>(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month  ,:,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(month+1,:,:,:)) ;     
  end
  elseif vec(1)==4
      if month==1
          var = 0.25*squeeze(var(4,:,:,:))+0.75*squeeze(var(1,:,:,:)) ;
      end
      if month==2
          var = squeeze(var(1,:,:,:)) ;
      end
      if month==3
          var = 0.75*squeeze(var(1,:,:,:))+0.25*squeeze(var(2,:,:,:)) ;
      end
      if month==4
          var = 0.25*squeeze(var(1,:,:,:))+0.75*squeeze(var(2,:,:,:)) ;
      end
      if month==5
          var = squeeze(var(2,:,:,:)) ;
      end
      if month==6
          var = 0.75*squeeze(var(2,:,:,:))+0.25*squeeze(var(3,:,:,:)) ;
      end
      if month==7
          var = 0.25*squeeze(var(2,:,:,:))+0.75*squeeze(var(3,:,:,:)) ;
      end
      if month==8
          var = squeeze(var(3,:,:,:)) ;
      end
      if month==9
          var = 0.75*squeeze(var(3,:,:,:))+0.25*squeeze(var(4,:,:,:)) ;
      end
      if month==10
          var = 0.25*squeeze(var(3,:,:,:))+0.75*squeeze(var(4,:,:,:)) ;
      end
      if month==11
          var = squeeze(var(4,:,:,:)) ;
      end
      if month==12
          var = 0.75*squeeze(var(4,:,:,:))+0.25*squeeze(var(1,:,:,:)) ;
      end
  end
  var = permute([permute(var(:,:,181:end),[1,3,2]) permute(var(:,:,1:180),[1,3,2])],[1,3,2]);
  var = var(:,j0:j1,i0:i1);
 
  disp('WARNING : fix for woa : bottom')
  for i=1:i1-i0+1
  for j=1:j1-j0+1
      indnan = min(find(isnan(squeeze(var(:,j,i))))) ;
      if (indnan~=1)
      var(indnan:nz,j,i) = var(indnan-1,j,i);
      end
  end
  end
  var(isnan(var))=0;  
  
  
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


   
