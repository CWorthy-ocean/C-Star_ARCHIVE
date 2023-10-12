%---------------------------------------------------------------------------------------
%
%  make_s2r
%
%  Generate boundary perimeter file from WOA and SSH  data.
%
%  Note that when run this script it tests for the presence of a .mat file
%  which contains various interpolation coefficients related to your child
%  and parent grids.  If the .mat file is not there it will calculate the coefficients
%
%
%  Jeroen Molemaker and Evan Mason in 2007-2009 at UCLA
%
%---------------------------------------------------------------------------------------
clear all
close all
disp(' ')

day_start = 1;
day_end = 31;

for yyyy=2012:2012

    clear list_soda_hbl
    clear list_soda_tr
    clear list_soda_vel
    clear list_soda_ssh
    clear BRYtime

%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS START HERE
%---------------------------------------------------------------------------------------
%
%---------------------------------------------------------------------------------------
% 1.  GENERAL
%---------------------------------------------------------------------------------------

    %  Data climatology file names:
    glorys_dir = ['/glade/scratch/bachman/GLORYS/NA/']
    glorys_file = '/mercatorglorys12v1_gl12_mean_'
    glorys_mon_tr = [glorys_dir, num2str(yyyy), glorys_file num2str(yyyy) '*.nc'];
    glorys_mon_vel= [glorys_dir, num2str(yyyy), glorys_file num2str(yyyy) '*.nc'];
    glorys_mon_ssh= [glorys_dir, num2str(yyyy), glorys_file num2str(yyyy) '*.nc'];


    %%%%%
    % ROMS info

    grdname = ['/glade/scratch/bachman/ROMS_tools/setup_s2r_phys_only/1.Make_grid/Wales0_grd.nc'];
    bryname = ['/glade/scratch/bachman/ROMS_tools/setup_s2r_phys_only/1.Make_grid/Wales0_bry_' num2str(yyyy) '.nc'];

    pars.theta_s = 5.0;
    pars.theta_b = 2.0;
    pars.hc     = 300.0;
    pars.N      = 100;
    pars.scoord = 'new2012';    % child 'new' or 'old' type scoord

    dateref = datenum(2000,1,1)

    % OBC flags and time
    obcflag        = [1 1 1 1];    % open boundaries flag (1=open , [S E N W])
    bry_cycle      = 0 %365.25;       % 0 means no cycle
    %mdays = [31 28 31 30 31 30 31 31 30 31 30 31] ;
    DT             = 1 ;           % time step of the forcing file (days)
    extraband      = 0  ;          % extra value at the boundary
%
%---------------------------------------------------------------------------------------
% USER-DEFINED VARIABLES & OPTIONS END HERE
%---------------------------------------------------------------------------------------
%
    rep_tr  = dir(glorys_mon_tr ) ;
    rep_vel = dir(glorys_mon_vel) ;
    rep_ssh = dir(glorys_mon_ssh) ;

    rep_tr = rep_tr(day_start:day_end);
    rep_vel = rep_vel(day_start:day_end);
    rep_ssh = rep_ssh(day_start:day_end);

    for l=1:length(rep_tr)
      list_soda_tr {l} = [glorys_dir num2str(yyyy) '/' rep_tr(l,1).name ];
      list_soda_vel{l} = [glorys_dir num2str(yyyy) '/' rep_vel(l,1).name];
      list_soda_ssh{l} = [glorys_dir num2str(yyyy) '/' rep_ssh(l,1).name];
    end

    if (extraband>0)
        if extraband>1
            disp('ERROR extraband>1 not coded') ; error ;
        end
        file_prev = dir([glorys_dir, num2str(yyyy-1) glorys_file num2str(yyyy-1) '*.nc']) ;
        file_prev = [glorys_dir, num2str(yyyy-1) '/'  file_prev(end,1).name];
        file_next = dir([glorys_dir, num2str(yyyy+1) glorys_file num2str(yyyy+1) '*.nc']) ;
        file_next = [glorys_dir, num2str(yyyy+1) '/' file_next(1,1).name];
        list_soda_ssh = [file_prev list_soda_ssh file_next];

        file_prev = dir([glorys_dir, num2str(yyyy-1) glorys_file num2str(yyyy-1) '*.nc']) ;
        file_prev = [glorys_dir, num2str(yyyy-1) '/' file_prev(end,1).name];
        file_next = dir([glorys_dir, num2str(yyyy+1) glorys_file, num2str(yyyy+1) '*.nc']) ;
        file_next = [glorys_dir, num2str(yyyy+1) '/' file_next(1,1).name];
        list_soda_tr = [file_prev list_soda_tr file_next];

        file_prev = dir([glorys_dir, num2str(yyyy-1) glorys_file num2str(yyyy-1) '*.nc']) ;
        file_prev = [glorys_dir, num2str(yyyy-1) '/' file_prev(end,1).name];
        file_next = dir([glorys_dir, num2str(yyyy+1) glorys_file num2str(yyyy+1) '*.nc']) ;
        file_next = [glorys_dir, num2str(yyyy+1) '/' file_next(1,1).name];
        list_soda_vel = [file_prev list_soda_vel file_next];
    end


    %%% BRY time
    BRYtime.mdays = [31 29 31 30 31 30 31 31 30 31 30 31] ;

    for t=1:length(list_soda_tr)
        file1 = list_soda_tr{t};

        filedate = regexp(file1,'\d\d\d\d\d\d\d\d','Match');
        filedate = filedate{1};

        BRYtime.year(t) = str2num(filedate(1:4));
        BRYtime.month(t) = str2num(filedate(5:6));
        BRYtime.day(t) = str2num(filedate(7:8));

        BRYtime.time(t) = datenum(BRYtime.year(t), BRYtime.month(t), BRYtime.day(t)) - dateref;

    end

    BRYtime.tstart = BRYtime.time(extraband+1);
    BRYtime.tend   = BRYtime.time(length(BRYtime.time)-extraband);
    BRYtime.cycle = bry_cycle ;


    % Create the bry file
    if ~exist(bryname)
      disp(['Creating boundary file: bryname' bryname]);
      r2r_create_bry(bryname,grdname,obcflag,pars,BRYtime);
    end

    %for days = 1:length(list_soda_tr)
        %disp(['%%----> Working on days ' num2str(days) '/' num2str(length(list_soda_tr)) ' <----%%'])
        %s2r_hv_glorys2(list_soda_tr{days},list_soda_vel{days},list_soda_ssh{days},grdname,bryname,days,pars,obcflag);
    %end
    s2r_hv_glorys(list_soda_tr,list_soda_vel,list_soda_ssh,grdname,bryname,days,pars,obcflag);

end
