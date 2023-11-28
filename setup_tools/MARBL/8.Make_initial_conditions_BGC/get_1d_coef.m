function [coef1d,elem1d] = get_1d_coef(zp,zc)

%
% Inputs:
%          parent z (zp)
%          child  z (zc)
%
% Ouputs:
%          elem - pointers to 1d gridded data (at zp locations) from
%                 which the interpolation is computed (1-2 for each child point)
%          coef - linear interpolation coefficients
%
%  Use:
%          To subsequently interpolate data from Fp(lonp,latp) to Fc(lonc,latc), the following
%          will work:      Fc  = sum(coef.*Fp(elem),3);  This line  should come in place of all
%          griddata calls. Since it avoids repeated triangulations and tsearches (that are done
%          with every call to griddata) it should be much faster.
%
%   Jeroen Molemaker, UCLA  2007


     [Np dum] = size(zp);
     [Nc dum] = size(zc);

     coef1d = zeros(Nc,2);
     elem1d = ones(Nc,2);

     ip = 1;
     for ic = 1:Nc
       while (zp(ip)<zc(ic)) & (ip<Np)
        ip = ip+1;
       end
       if ip == 1
         coef1d(ic,1) = 1.0;
         elem1d(ic,1) = ip;
	 continue
       end
       if zp(ip) < zc(ic)
         coef1d(ic,1) = 1.0;
         elem1d(ic,1) = ip;
	 continue
       end
       alp = (zc(ic)-zp(ip-1))/(zp(ip)-zp(ip-1));
       coef1d(ic,1) = alp;
       elem1d(ic,1) = ip;
       coef1d(ic,2) = 1-alp;
       elem1d(ic,2) = ip-1;
     end
