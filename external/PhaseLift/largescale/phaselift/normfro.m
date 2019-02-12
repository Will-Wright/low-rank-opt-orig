function nf=normfro(A,B,r)
% NORMFRO  - Frobenius norm of matrix A (or of A-B)
%
% Usage:     nf = normfro(A,B,1);
%            if flag 1 is used, norm(A-B,'fro')/norm(A,'fro') is computed
%
%
% Input:     A,B matrices of same size
% Output:     
%
% See also:   

% Author:     T.Strohmer 02.Jan.1996 (strohmer@tyche.mat.univie.ac.at)
% Modified:   
%
% Literature: 
% 
% Copyright:  (c) NUHAG, Dept. Math., University of Vienna, Austria
%             Permission is granted to modify and re-distribute this
%             code in any manner as long as this notice is preserved.
%             All standard disclaimers apply.
% 
% Externals:  


if nargin == 1
  nf = norm(A,'fro');
elseif nargin == 2
  nf = norm(A-B,'fro');
elseif nargin == 3
  nf = norm(A-B,'fro')/norm(A,'fro');
end

% end of file
