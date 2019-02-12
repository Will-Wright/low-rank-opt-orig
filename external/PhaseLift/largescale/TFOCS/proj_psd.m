function op = proj_psd( LARGESCALE, K )

% PROJ_PSD  Positive semidefinite cone.
%   OP = PROJ_PSD() returns a function that implements
%   the projection onto the semidefinite cone:
%	X = argmin_{min(eig(X))>=0} norm(X-Y,'fro')
%
%   OP = PROJ_PSD( LARGESCALE )
%     performs the same computation, but in a more efficient
%     manner for the case of sparse (and low-rank) matrices
%
%   OP = PROJ_PSD( LARGESCALE, k )
%     only returns at most a rank k matrix

if nargin == 0, LARGESCALE = false; end
if nargin < 2,  K = Inf;  end

if ~LARGESCALE
    op = @proj_psd_impl;
else
    op = @(varargin)proj_psd_largescale( K, varargin{:} );
end


function [ v, X ] = proj_psd_impl( X, t )
if nargin > 1 && t > 0,
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


function [ v, X ] = proj_psd_largescale(K,X, t )
if nnargin > 2 && t > 0,
	v = 0;
    [V,D]=eig(X+X');
    D  = max(0.5*diag(D),0);
    tt = D > 0;
    V  = bsxfun(@times,V(:,tt),sqrt(D(tt,:))');
    X  = V * V';
else
%     s = eig(X);
%     if min(s) < -8*eps*max(s),
%         v = Inf;
%     else
%     	v = 0;
%    	end

    % This is potentially expensive.
    
    
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

