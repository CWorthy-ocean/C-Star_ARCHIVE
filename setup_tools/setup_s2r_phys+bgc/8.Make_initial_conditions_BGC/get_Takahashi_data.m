function [var,loni,lati] = get_Takahashi_data(Takahashi_file,variable,month,lon,lat,BRYtime,days);

% Processes and returns data from woa18 monthly climatology that is mask filled
% and extended from -180 to 360 degrees.

%%%% FIND the file to be read %%%%%%
  infile = Takahashi_file ;
  disp(['get_Takahashi : reading file : ' infile])

%%%% read lon/lat %%%%%%
  loni = ncread(infile,'lon')';
  lati = ncread(infile,'lat')';
  %%%i0 = find(loni<min(min(lon)),1,'last');
  %%%i1 = find(loni>max(max(lon)),1,'first');
  %%%j0 = find(lati<min(min(lat)),1,'last');
  %%%j1 = find(lati>max(max(lat)),1,'first');
  [loni,lati] = meshgrid(loni,lati);
 %figure(1);plot(lon,lat,'.k')
 %figure(2);plot(loni,lati,'.k')
% error

  %%%loni = loni(j0:j1,i0:i1);
  %%%lati = lati(j0:j1,i0:i1);

  [ny,nx] = size(loni);

  if strcmp(variable,'pCO2')==1
      name='pco2_sw';
  else
      name=variable;
  end

  var = permute(squeeze(ncread(infile,name)),[3 2 1]);
%   var = squeeze(var(month,:,:)) ;
  if ( (month == 1) && BRYtime.day(days)<=(BRYtime.mdays(month)/2) )
      var =  BRYtime.interp_ratio(days) *squeeze(var(month,:,:)) + ...
          (1-BRYtime.interp_ratio(days))*squeeze(var(   12,:,:)) ;
  elseif ( (month == 1) && BRYtime.day(days)>(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(    2,:,:)) ;
  elseif ( (month == 12) && BRYtime.day(days)<=(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(   11,:,:)) ;
  elseif ( (month == 12) && BRYtime.day(days)>(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(    1,:,:)) ;
  elseif ( BRYtime.day(days)<=(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month  ,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(month-1,:,:)) ;
  elseif ( BRYtime.day(days)>(BRYtime.mdays(month)/2) )
       var =  BRYtime.interp_ratio(days) *squeeze(var(month  ,:,:)) + ...
           (1-BRYtime.interp_ratio(days))*squeeze(var(month+1,:,:)) ;
  end

  %%%var = var(j0:j1,i0:i1);


end



