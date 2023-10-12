function alk  = lee_sfc_alk(sst, sss, lon, lat, basin)
%LEE_SFC_ALK Calculate surface Alkalinity from SST and SSS
%   USAGE: alk = lee_sfc_alk(sst, sss, lon, lat, basin)
% After Lee et al.: Global relationships of total alkalinity with salinity and temperature
% in surface waters of the world's oceans, GRL, VOL. 33, L19605, doi:10.1029/2006GL027207, 2006)
%
% Reimplemented by M. Munnich 7/2011

% Departures from reference salinity and temperature
dS = sss-35.0;
dS2 = dS.^2;
dT = sst-20.0;
dT2 = dT.^2;
lon = mod(lon,360);
%
% Zone 1 (Sub)tropics
%
alk1 = 2305 + 58.66*dS + 2.32*dS2 -  1.41*dT + 0.040*dT2;
%
% Zone 2 (Equatorial Upwelling Pacific)
%
alk2 = 2294 + 64.88*dS + 0.39*dS2 - 4.52*(sst - 29) - 0.232*(sst - 29).^2;
%
% Zone 3 North Atlantic
%
alk3 = 2305 + 53.97*dS + 2.74*dS2 - 1.16*dT - 0.040*dT2;
%
% Zone 4 North Pacific 
%
alk4 = 2305 + 53.23*dS + 1.85*dS2 - 14.72*dT - 0.158*dT2 + 0.062*dT.*lon;
%
% Zone 5 Southern Ocean ("used as default below")
%
alk5 = 2305 + 52.48*dS  + 2.85*dS2 - 0.49*dT + 0.086*dT2;

alk = alk5;

%
% Overwrite Souther Ocean default based on lat, lon, basin and temperature
%
% Subtropics
%
subtrop = sst>=20;
alk(subtrop) = alk1(subtrop);
%
% North Pacific
%
%% North of the Equator the North Pacific is adjacent to the Tropics in the Pacific
northpac = lat > 0 & sst <= 20 & basin == 2 ; 
alk(northpac) = alk4(northpac);
%
% North Atlantic/Mediterrian (basin == 4)
%
%% North of the Equator the North Atlantic is adjacent to the Tropics in the Atlantic
northatl = lat > 0 & sst <= 20 & ( basin == 1 | basin == 4 );
alk(northatl) = alk3(northatl);
%
% Equatorial Pacific
%
eqpac = ( abs(lat)<20 | abs(lat)<10+(lon-220)/3 & lon<250) ...
             & lon>=220 & sst<=29 & basin==2;
%eqpac = (  lon> 220 & abs(lat)<20+(lon-220)/6  ... % lon=250 -> lat<25
%        ... % One has to decide which is the neighboring region of the
%        ... % TropicsRegion EquatorPacific or NorthPacifi/SouthernOcean 
%        ... % For SST < 29 in the TropicsRegion.
%        ... % Here we interpret it in spirit of what is close to EquatorialRegion
%        ... % upwelling NOT exactly what is closer in distance.
%          | lon > 180 & abs(lat)< (lon-180)/2  ) ...  % lon = 220 -> lat<20 
%        & sst<=29 & basin==2 & abs(lat)<25;
alk(eqpac) = alk2(eqpac);
% Correct low temperature tropical values near eq upwelling region
%eqpac = lon>220 & abs(lat)<25 & basin==2 & abs(lat)<10

%eqtrop = 
%
% Smoothing between across transition lines  between regions
% Tsmooth: half width of smoothing transition in temperature
%
% e.g.: Tsmooth = 1:  20C transition is linearily smoothed between T = 19C
% and T = 21C
%
Tsmooth = 1;
if Tsmooth > 0
    [NX,NY] = size(sst);
    % 20C, 29C: Lee et al transition temperatures
    T1 = 20-Tsmooth;
    T2 = 20+Tsmooth;
    T3 = 29-Tsmooth;
    T4 = 29+Tsmooth;
    dTi = 1/(2*Tsmooth);
    for i=1:NX
        for j = 1:NY
            if eqpac(i,j)
                if sst(i,j) > T3 && sst(i,j)<T4
                    alk(i,j) = alk1(i,j)*dTi*(sst(i,j)-T3)+alk2(i,j)*dTi*(T4-sst(i,j));
                end
            elseif lat(i,j) <=0 && sst(i,j)>=T1 && sst(i,j)<=T2 % Southern Hemisphere
                alk(i,j) = alk5(i,j)*dTi*(T2-sst(i,j))+alk1(i,j)*dTi*(sst(i,j)-T1);
            elseif lat(i,j)>=0 && sst(i,j)>=T1 && sst(i,j) <= T2 && basin(i,j) == 2 % Northern Hemisphere Pacific
                alk(i,j) = alk4(i,j)*dTi*(T2-sst(i,j))+alk1(i,j)*dTi*(sst(i,j)-T1);
            elseif lat(i,j)>=0 && sst(i,j)>=T1 && sst(i,j) <= T2 && (basin(i,j) == 1 || basin(i,j) == 4) % Northern Hemisphere Atlantic
                alk(i,j) = alk4(i,j)*dTi*(T2-sst(i,j))+alk1(i,j)*dTi*(sst(i,j)-T1);
            end
        end
    end
end
%
% Smoothing latitudinally of eqPacific/Tropic boundary
latsmooth = 2; 
if latsmooth > 0
    if any(any(basin == 2 & abs(lat)<20+latsmooth))
        for i=1:NX
            for j = 1:NY
                if basin(i,j) == 2
                if abs(lat(i,j)) < 20+latsmooth && lon(i,j) > 220 
                    if lon(i,j) > 250 
                        latbry = 20;
                    elseif lon(i,j)<=250
                        latbry = 10+(lon(i,j)-220)/3;
                    end
                    dlat = 0.5*(abs(lat(i,j))-(latbry-latsmooth))/latsmooth;
                    if dlat <=1 && dlat > 0
                        alk(i,j) = alk1(i,j)*dlat+alk2(i,j)*(1-dlat);
                    end
                end
                end
            end
        end
    end
end
                

end



