function fm = fillmask(f,type,maskp,nnel);

%     Put values in the masked points that are based on nearest
%     non-masked neighbour or zero (type = 1/0)
%     This take the place of 'fill missing value' and 'nonan'
%
%     input:  f      2d or 3d grid of values to be mask filled
%             type   mask extension type: 0 zeros/1 nearest neighbor
%                    extrapolation
%     output: fm     2d or 3d grid of values that are filled
%
%     (c) 2007, Jeroen Molemaker,  UCLA

    if ndims(f)==3
      % 3d grid
      [n m l] = size(f);

      if type==1   %% nearest neighbor extrapolate
        for k = 1:n
          fk  = f(k,:,:);
          fk(find(maskp==0)) = fk(nnel(find(maskp==0)));
          fm(k,:,:) = fk;
        end
      else          %% zeros
        for k = 1:n
          fk = f(k,:,:);
          fk(find(maskp==0)) = 0.0*fk(find(maskp==0));
          fm(k,:,:) = fk;
        end
      end

    else
      % 2d grid
      if type==1   %% nearest neighbor extrapolate
          f(find(maskp==0)) = f(nnel(find(maskp==0))); 
          fm = f;
      else          %% zeros
          f(find(maskp==0)) = 0.0*f(find(maskp==0));
          fm = f;
      end
    end

    return
