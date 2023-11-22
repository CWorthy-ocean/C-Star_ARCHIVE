function s2r_hvbgc_recalc(BGC_INI,grdname,bryname,chdscd,obcflag);
%--------------------------------------------------------------
%  Produce a ROMS boundary file from SODA 2.0.4 data
%
%  Inspired by Roms_tools (IRD).
%  Thanks to Pierrick, Patrick (IRD), Francois (UCLA), Yusuke (UCLA)
%  Jeroen Molemaker (UCLA); nmolem@ucla.edu
%--------------------------------------------------------------

% Get S-coordinate params for child grid
  theta_b_c = chdscd.theta_b;
  theta_s_c = chdscd.theta_s;
  hc_c      = chdscd.hc;
  N_c       = chdscd.N;
  scoord_c  = chdscd.scoord;


  [npc mpc] = size(ncread(grdname,'h'));


 for bnd = 1:4

  disp('-------------------------------------------------------------')
  if ~obcflag(bnd)
    disp('Closed boundary')
    continue
  end
  if bnd==1
   disp('South boundary')
   i0 = 1;
   i1 = npc;
   j0 = 1;
   j1 = 2;
   suffix='_south';
  end
  if bnd==2
   disp('East boundary')
   i0 = npc-1;
   i1 = npc;
   j0 = 1;
   j1 = mpc;
   suffix='_east';
  end
  if bnd==3
   disp('North boundary')
   i0 = 1;
   i1 = npc;
   j0 = mpc-1;
   j1 = mpc;
   suffix='_north';
  end
  if bnd==4
   disp('West boundary')
   i0 = 1;
   i1 = 2;
   j0 = 1;
   j1 = mpc;
   suffix='_west';
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[C,trc]=find(ismember(BGC_INI.bgc_tracer,'DIC')==1) ;
%   if (C==1)
%       disp('Computing DIC and Alk from surface pCO2')
%       [srf_NTA_vol,srf_dic_vol]=det_sfc_dic_nta_BRY(bryname,grdname,chdscd,bnd) ;
%       dic_alk_profile_BRY(bryname,grdname,chdscd,srf_NTA_vol,srf_dic_vol,-750,bnd) ;
%   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   vec = strfind(BGC_INI.bgc_frctype,'s') ;
   for trc=1:length(vec)
       if vec{trc}==1
       scaled_tracer = BGC_INI.bgc_tracer{trc} ;
       ref_ind       = str2num(BGC_INI.bgc_frctype{trc}(3:end)) ;
       tracer_scaled = BGC_INI.bgc_tracer{ref_ind} ;
       factor        = str2num(BGC_INI.bgc_frcini{trc}) ;
       disp(['SCALING variable ' scaled_tracer ' from ' tracer_scaled ' and factor ' num2str(factor)])
       var = ncread(bryname,[tracer_scaled suffix]) ;
       ncwrite(bryname,[scaled_tracer suffix],var*factor,[1 1 1]);
       end
   end


end    % End loop bnd


  return



















