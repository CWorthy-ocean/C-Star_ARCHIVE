
frcname = '/glade/cheyenne/scratch/bachman/C-Star/setup_tools/MARBL/9.Make_boundary_forcing_BGC/roms_bryBGC.nc';

disp(['Changing bry_time in ' frcname])

time = ncread( frcname , 'bry_time' ) ;

time = 4383 + [0 31 60 91 121 152 182 213 244 274 305 335] ;
ncwrite(frcname , 'bry_time' , time ) ;

