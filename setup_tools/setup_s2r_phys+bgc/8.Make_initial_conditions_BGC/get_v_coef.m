function Av = get_v_coef(zp,zc);

    [Np Mc Lc] = size(zp);
    [Nc Mc Lc] = size(zc);

    ca = zeros(Nc,Mc,Lc,2);
    ic = zeros(Nc,Mc,Lc,2);
    jc = zeros(Nc,Mc,Lc,2);
    for i = 1:Mc
     for j = 1:Lc
      ka = [sub2ind([Nc Mc Lc],1,i,j):sub2ind([Nc Mc Lc],Nc,i,j)];
      k2d = repmat(ka,2,1)';
      [coef1d,elem1d] = get_1d_coef(squeeze(zp(:,i,j)),squeeze(zc(:,i,j)));
      ca(:,i,j,:) = coef1d;
      ic(:,i,j,:) = k2d;
      jc(:,i,j,:) = elem1d + Np*(i-1) + Np*Mc*(j-1);
      end
    end
    ndimc = Nc*Mc*Lc;
    ndimp = Np*Mc*Lc;
    ic = reshape(ic,ndimc*2,1);
    jc = reshape(jc,ndimc*2,1);
    ca = reshape(ca,ndimc*2,1);
    Av = sparse(ic,jc,ca,ndimc,ndimp);

    return


    ic1 = [1:Nc];
    ic2 = reshape(repmat(ic1,2,1),2*Nc,1);
    jc2 = reshape(elem1d',2*Nc,1);
    ca2 = reshape(coef1d',2*Nc,1);
    A = sparse(ic2,jc2,ca2,Nc,Np);

    fc = sum(coef1d.*fp(elem1d),2);
    fc0= A*fp';
    fc1= interp1(zp,fp,zc);
    fc2= interp1(zp,fp,zc,'spline');
    plot(zp,fp,'o')
    hold on
    plot(zc,fc1,'g')
    plot(zc,fc)
    plot(zc,fc2,'r')
    hold off

