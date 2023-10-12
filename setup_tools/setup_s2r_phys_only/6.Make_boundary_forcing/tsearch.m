function k = tsearch(x,y,tri,xi,yi)
%TSEARCH Closest triangle search.
%   T = TSEARCH(X,Y,TRI,XI,YI) returns the index of the enclosing Delaunay
%   triangle for each point in XI,YI so that the enclosing triangle for
%   point (XI(k),YI(k)) is TRI(T(k),:).  TSEARCH returns NaN for all
%   points outside the convex hull.  Requires a triangulation TRI of the
%   points X,Y obtained from DELAUNAY.
%
%   See also DelaunayTri, DELAUNAY, DSEARCH, QHULL, TSEARCHN, DELAUNAYN.

%   Relies on the MEX file tsrchmx to do most of the work.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 1.17.4.4 $  $Date: 2008/12/01 07:18:40 $


    % Create TriangleGraph which is a sparse connectivity matrix
    ntri = size(tri,1);
    triangle_graph = sparse( repmat((1:ntri)',1,3), tri, 1, ntri, numel(x));
    
    % call tsrchmx to do the work
    %X, Y are treated as vectors inside mex function
    k = tsrchmx(x,y,tri,xi,yi,triangle_graph);


