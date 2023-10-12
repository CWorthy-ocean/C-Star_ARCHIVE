function chlo=extr_chlo(Cpd,z,convert)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function chlo=extr_chlo(Cpd,z);
%
%  pierrick 2001
%
%  vertical extrapolation of chlorophyll from the surface values
%  using Morel and Berthon (1989) parameterization
%  ref:  Morel and Berthon, Surface pigments, algal biomass
%        profiles, and potential production of the euphotic layer:
%        Relationships reinvestigated in view of remote-sensing 
%        applications. Limnol. Oceanogr., 34, 1989, 1545-1562.
%
%    input:
%
%  Cpd : mean pigment concentration within the surface layer
%        2D horizontal matrix. UNITS: (mg Chla/m3) if convert=0
%        (default), otherwise (micro mole/l)
%  z   : vertical positions (m, negative) 3D matrix
%  convert (OPTIONAL argument) : convert from umol/l to mgChla/m3 if 1
%                                DO NOT CHANGE UNITS IF 0, i.e. if
%                                Cpd is already given in mgChla/m3
%
%    output:
%
%  chlo: chlorophyll concentration (mgC.m-3) 3D matrix (convert = 1)
%  chlo: chlorophyll concentration (mgChla.m-3) 3D matrix (convert = 0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global VERBOSE
if isempty(VERBOSE)
    VERBOSE = 0;
end
if VERBOSE > 2
    disp('Morel Berthon vertical chlorophyll distribution.')
end
% Note: We keep depth as the 1st (fasted changing dimension)
%  
sz = size(z);

Cpd = repmat(reshape(Cpd,[1,sz(2:end)]),[sz(1),1,1]);
%
% total pigment content within the euphotic layer (Ctot)
% expressed as the mean pigment concentration within 
% the surface layer (Cpd) (stratified waters)
%
Ctot=38*Cpd.^0.425; % Eq. (2b) in Morel & Berthon (1989)
Ctot(Cpd>1)=40.2*Cpd(Cpd>1).^0.507; % Eq. (2c) in Morel & Berthon (1989)
%
% depth of the euphotic layer
%
Ze=-568.2*Ctot.^(-0.746); % Eq. (1a) in Morel & Berthon (1989)
Ze(Ze<=-102)=-200*Ctot(Ze<=-102).^(-0.293); % Eq. (1b) in Morel & Berthon (1989)
%mm zeta=z./Ze; %mm not needed to find Ze
%
% mean pigment concentration within 
% the euphotic layer (Cze)
% espressed as the mean pigment concentration within 
% the surface layer (Cpd) (stratified waters)
%
Cze=1.12*Cpd.^0.803; % Eq. on top left of p. 1557 in Morel & Berthon (1989)
%
% vertical chlorophyll profile 
%
lc=log10(Cpd);
clear Cpd
lc2=lc.^2;
lc3=lc2.*lc;
% equations after (6) on p. 1557 in Morel & Berthon (1989)
Cb=0.768+0.087*lc-0.179*lc2-0.025*lc3;
Cmax=0.299-0.289*lc+0.579*lc2;
zetamax=0.600-0.640*lc+0.021*lc2+0.115*lc3;
dzeta=0.710+0.159*lc+0.021*lc2;
clear lc lc2 lc3
%dzetai = 1.0./dzeta;
%
% recompute iteratively the euphotic layer depth
%
c1=-(Cze.*Cb+Cze.*Cmax.*0.5.*sqrt(pi).*dzeta.*...
           (erf((1-zetamax)./dzeta)-erf(-zetamax./dzeta)));
eps=1000;
while eps>0.01
  Zeold=Ze;
  %mm Ctot=-Ze.*(Cze.*Cb+Cze.*Cmax.*0.5.*sqrt(pi).*dzeta.*...
  %mm           (erf((1-zetamax)./dzeta)-erf(-zetamax./dzeta))); 
  %mm above is Morel Berthon (1989) eq (6) integrated over z. 
  %mm now rewritten for speed-up using c1
  Ctot=Ze.*c1;
  Ze=-568.2*Ctot.^(-0.746); % Eq. (1a) in Morel & Berthon (1989)
  Ze(Ze<=-102)=-200*Ctot(Ze<=-102).^(-0.293); % Eq. (1b) in Morel & Berthon (1989)
  Ze=0.5*(Ze+Zeold);
  eps=max(max(abs(Zeold-Ze)));
end
clear c1 Ctot 
  zeta=z./Ze; %mm put outside while loop
% the second line is eq. (6) in Morel & Berthon (1989)
chlo=0.5.*(1+tanh(2.*(2.*Ze-z)./Ze)).*...
  Cze.*(Cb+Cmax.*exp(-((zeta-zetamax).*dzeta).^2));
chlo(chlo<0)=0;
end
%
