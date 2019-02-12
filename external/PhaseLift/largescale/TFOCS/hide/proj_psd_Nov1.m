function [ v, X ] = proj_psd( X, t )
% PROJ_PSD  Positive semidefinite cone.
%   OP = PROJ_PDS() returns a function that implements
%   the projection onto the semidefinite cone:
%	X = argmin_{min(eig(X))>=0} norm(X-Y,'fro')
%
% This uses a naive approach that will be more expensive than necessary 
% for large matrices that are known to be low rank. A future implementation
% of TFOCS will be able to handle low-rank matrices more efficiently.

if nargin == 0,
	v = @proj_psd;
elseif nargin > 1 && t > 0,
	v = 0;
    [V,D]=eig(X+X');
    D  = max(0.5*diag(D),0);
    tt = D > 0;
    V  = bsxfun(@times,V(:,tt),sqrt(D(tt,:))');
    X  = V * V';
else
    s = eig(X);
    if min(s) < -8*eps*max(s),
        v = Inf;
    else
    	v = 0;
   	end
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

