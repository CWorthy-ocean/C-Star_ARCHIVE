gridfile = '/glade/scratch/bachman/ROMS_tools/Iceland0_BGC2/1.Make_grid/Iceland0_grd.nc';

mask = ncread(gridfile,'mask_rho');
%mask = ncread(gridfile,'mask');
[ny,nx] = size(mask);

reg = bwlabel(mask,4);

lint =  0; %% size of largest region
lreg = 0; %% number of largest region
nreg = max(max(reg)); %% number of regions
for i = 1:nreg
  int = sum(sum(reg==i));
  if int>lint
    lreg = i;
    lint = int;
  end
end


for ireg = 1:nreg
  if ireg~=lreg
    %% check before setting to zero
    int = sum(sum(reg==ireg));
    if int>nx*ny/10
      disp(['region: ' num2str(ireg) 'is large.'])
    else
     mask(reg==ireg) = 0;
    end
  end
end

if 1
  subplot(1,2,1)
  imagesc(reg);axis xy;colorbar
  subplot(1,2,2)
  imagesc(mask);axis xy;colorbar
end
  ncwrite(gridfile,'mask_rho',mask);


