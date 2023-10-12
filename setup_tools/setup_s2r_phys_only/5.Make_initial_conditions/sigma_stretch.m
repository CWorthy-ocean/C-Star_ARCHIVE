function [sc,Cs] = sigma_stretch(theta_s,theta_b,N,type,sigmatype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGMA_STRECH: Compute S-coordinate stretching factor Cs
% 
%   USAGE:
%   [sc,Cs] = sigma_stretch(theta_s,theta_b,N,type);
%
%  Input:
%    theta_s, theta_b, hc: stretching parameters
%    N: # depth layers
%    type: grid type   'r': rho point 'w': w point 
%
%  Output:
%
%    Cs   S-coordinate stretching factor Cs 
%    sc   level location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global SIGMA_COORD_TYPE
if nargin == 5
    SIGMA_COORD_TYPE = sigmatype;
end
%
% Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
cff1=1./sinh(theta_s);
cff2=0.5/tanh(0.5*theta_s);
if type=='w'
  sc=((0:N)-N)/N;
  % N=N+1;
else
  sc=((1:N)-N-0.5)/N;
end
if (SIGMA_COORD_TYPE <= 2)
	Cs=(1.-theta_b)*cff1*sinh(theta_s*sc)...
    	+theta_b*(cff2*tanh(theta_s*(sc+0.5))-0.5);
elseif (SIGMA_COORD_TYPE == 3)
	Cs=cs_f(sc,theta_s,theta_b);
else
	error('Unknown Sigma coord type')
end
return

function cs=cs_f(sc, theta_s, theta_b)
if (theta_s > 0.)
    csrf=(1.-cosh(theta_s*sc))/(cosh(theta_s)-1.);
else
    csrf=-sc.^2;
end
if (theta_b > 0.)
    cs=(exp(theta_b*csrf)-1.)/(1.-exp(-theta_b));
else
    cs=csrf;
end
return

