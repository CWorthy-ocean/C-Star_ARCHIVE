 function sms = fix_bry_mask(ename,obj_name,mskp)
 % Finds points that are in the parent mask and
 % but not masked in the child
 % and replace parent indices with nearest neighbor point

% disp(['Fixing: ', obj_name]);

  ijp = ncread(ename,obj_name);
  iv = ijp(:,1);
  jv = ijp(:,2);

  if strcmp(obj_name(end),'r')
    ivs = iv + 1.5;
    jvs = jv + 1.5;
    sms = iv*0;
    for i = 1:length(iv)
      ill = floor(ivs(i));
      jll = floor(jvs(i));
      sms(i) = sum(sum(mskp(ill:ill+1,jll:jll+1)));
    end   
  end

  if strcmp(obj_name(end),'u')|| strcmp(obj_name(end),'v')
    % due to rotations, we need to read both u and v for each.
    % So both of them need to be not fully in the mask. 
    ivu = iv + 1.0;
    jvu = jv + 1.5;
    msku = mskp(1:end-1,:).*mskp(2:end,:);
    smsu = iv*0;
    for i = 1:length(iv)
      ill = floor(ivu(i));
      jll = floor(jvu(i));
      smsu(i) = sum(sum(msku(ill:ill+1,jll:jll+1)));
    end

    ivv = iv + 1.5;
    jvv = jv + 1.0;
    mskv = mskp(:,1:end-1).*mskp(:,2:end);
    smsv = iv*0;
    for i = 1:length(iv)
      ill = floor(ivv(i));
      jll = floor(jvv(i));
      smsv(i) = sum(sum(mskv(ill:ill+1,jll:jll+1)));
    end   
    sms = smsu.*smsv;
  end
  
  x = 1:length(sms);
  ivm = iv;
  jvm = jv;
  xm = x;
  ivm(sms==0) = [];
  jvm(sms==0) = [];
  xm(sms==0) = [];

% indx = find((sms-1).*obj_msk<0);
  indx = find(sms==0);
  display(['Fixing ',num2str(length(indx)),' points that are inside parent mask'])
  iv(indx) = interp1(xm,ivm,x(indx),'nearest','extrap');
  jv(indx) = interp1(xm,jvm,x(indx),'nearest','extrap');

  ijp(:,1) = iv;
  ijp(:,2) = jv;

  ncwrite(ename,obj_name,ijp);
