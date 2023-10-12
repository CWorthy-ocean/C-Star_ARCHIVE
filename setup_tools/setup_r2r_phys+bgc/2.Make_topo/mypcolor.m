function mypcolor(x,y,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A shifted version of pcolor to get the positions correct.
%   Jeroen Molemaker UCLA 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [Mp,Lp]=size(x);
 M=Mp-1;
 L=Lp-1;

 x_pcol = zeros(Mp+1,Lp+1);
 y_pcol = zeros(Mp+1,Lp+1);
 f_pcol = zeros(Mp+1,Lp+1);
%
 x_tmp             = 0.5*(x(:,1:L)     + x(:,2:Lp)    );
 x_pcol(2:Mp,2:Lp) = 0.5*(x_tmp(1:M,:) + x_tmp(2:Mp,:));
 x_pcol(1,:)       = 2*x_pcol(2,:)     - x_pcol(3,:);
 x_pcol(:,1)       = 2*x_pcol(:,2)     - x_pcol(:,3);
 x_pcol(end,:)     = 2*x_pcol(end-1,:) - x_pcol(end-2,:);
 x_pcol(:,end)     = 2*x_pcol(:,end-1) - x_pcol(:,end-2);
 y_tmp             = 0.5*(y(:,1:L)     + y(:,2:Lp)    );
 y_pcol(2:Mp,2:Lp) = 0.5*(y_tmp(1:M,:) + y_tmp(2:Mp,:));
 y_pcol(1,:)       = 2*y_pcol(2,:)     - y_pcol(3,:);
 y_pcol(:,1)       = 2*y_pcol(:,2)     - y_pcol(:,3);
 y_pcol(end,:)     = 2*y_pcol(end-1,:) - y_pcol(end-2,:);
 y_pcol(:,end)     = 2*y_pcol(:,end-1) - y_pcol(:,end-2);

 f_pcol(1:Mp,1:Lp) = f;
 f_pcol(end,:)     = f_pcol(end-1,:);
 f_pcol(:,end)     = f_pcol(:,end-1);


 pcolor(x_pcol,y_pcol,f_pcol);shading flat;colorbar
return













