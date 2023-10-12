
function  [lonc,latc] = fix_outside_child(lonc,latc,t);

  [Mc,Lc]=size(lonc);
  disp('Fixing outside child points, make sure these are masked');

  %% move outside child points away:
  val_lonc = lonc;
  val_latc = latc;
  val_lonc(~isfinite(t)) = val_lonc(~isfinite(t)) + 1000;  %% move the outside points far away
  val_latc(~isfinite(t)) = val_latc(~isfinite(t)) + 1000;  %% move the outside points far away

  val_tri = delaunay(val_lonc,val_latc);
  Xc    = [reshape(    lonc,Mc*Lc,1) reshape(    latc,Mc*Lc,1) ];
  Xval  = [reshape(val_lonc,Mc*Lc,1) reshape(val_latc,Mc*Lc,1) ];
  nearel= dsearchn(Xval,val_tri,Xc);               %% find the nearest inside  neighbor
  val_lonc(~isfinite(t)) = val_lonc(nearel(~isfinite(t)));
  val_latc(~isfinite(t)) = val_latc(nearel(~isfinite(t)));
  lonc = val_lonc;
  latc = val_latc;

  return
