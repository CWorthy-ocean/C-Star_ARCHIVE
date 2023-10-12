

% ROMS parent and child grid directories
  %pdir = '/paracas/nmolem/PACHUG/';
  %pdir = '/glade/scratch/bachman/ROMS_tools/EASY/';
  cdir = '/glade/scratch/bachman/ROMS_tools/Iceland0_BGC2/1.Make_grid/';
  %pgrid = 'CT_parent.nc';
  GLORYS_bathy = '/glade/scratch/bachman/GLORYS/GLO-MFC_001_030_mask_bathy.nc';
  cgrid = 'Iceland0_grd.nc';

  %pgrid = [pdir pgrid]
  cgrid = [cdir cgrid]


   rmax = 0.2;
   hmin = 5;
   offset = 0.0;

   gridfile = cgrid;
   lsmooth

   %if exist('pgrid')
%  Match boundary topo: Only match to parent topography on open boundaries

%%% COMMENTING THIS OUT WHEN GLORYS IS THE PARENT

%    obcflag              = [1 1 1 1];      % open boundaries flag (1=open , [S E N W])
%    mod_cgrid
   %end






