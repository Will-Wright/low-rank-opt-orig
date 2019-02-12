function [err,c] = phaserr(x,xa)
% PHASERR - computes the l_2-error between two vectors or matrices 
%           modulo a phase factor 
%
% Usage:      err = phaserr(x,xa)
%
% Input:  
%             
% Output:    
%

% Author:     T.Strohmer (strohmer@math.ucdavis.edu)
%
% 
% Copyright:  (c) T.Strohmer, Dept. of Mathematics, University of 
%             California, Davis, USA.

ndx1 = find(abs(x(:)) > 1e-6);
ndx2 = find(abs(xa(:)) > 1e-6);
ndx = intersect(ndx1,ndx2);
c = mean(x(ndx)./xa(ndx));
c = c/abs(c);
xa = c*xa;
err = norm(x(:)-xa(:))/norm(x(:));





% end of file
