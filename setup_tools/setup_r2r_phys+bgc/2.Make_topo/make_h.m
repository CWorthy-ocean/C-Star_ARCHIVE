

% ROMS parent and child grid directories
  %pdir = '/paracas/nmolem/PACHUG/';
  pdir = '/glade/scratch/bachman/ROMS_tools/Iceland0/EASY/';
  cdir = '/glade/scratch/bachman/ROMS_tools/setup_r2r_phys+bgc/1.Make_grid/';
  pgrid = 'Iceland0_grd.nc';
  cgrid = 'Iceland1_grd.nc';

  pgrid = [pdir pgrid]
  cgrid = [cdir cgrid]


   rmax = 0.2;
   hmin = 5;
   offset = 0.0;

   gridfile = cgrid;
   lsmooth

   %if exist('pgrid')
%  Match boundary topo: Only match to parent topography on open boundaries

%%% COMMENTING THIS OUT WHEN GLORYS IS THE PARENT

    obcflag              = [1 1 1 1];      % open boundaries flag (1=open , [S E N W])
    mod_cgrid2
   %end






