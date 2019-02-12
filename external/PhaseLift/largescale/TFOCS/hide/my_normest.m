function [nrmX,cnt] = my_normest( varargin )
%MY_NORMEST Estimate the spectral norm.
%   MY_NORMEST(S) is an estimate of the spectral norm of the matrix S
%       using the power method.
%   [nrm,cnt] = MY_NORMEST(..) also gives the number of iterations used.
%
%   MY_NORMEST(S,St,n) uses the funcions S() and St() to compute
%   the forward and transpose multiplication (where S is m x n)
%   If St is the empty matrix, then we assume S = S'
%
%   MY_NORMEST(...,tol) uses relative error tol as stopping criteria
%       [default: 1e-6 ]
%   MY_NORMEST( ..., tol, maxiter ) uses at most maxiter iterations
%       [default: 20 ]


error(nargchk(1,5,nargin));

tol     = 1e-6;
maxiter = 20;
base    = 1;
S = varargin{1};
if isa(S,'function_handle')
    error(nargchk(2,5,nargin));
    IMPLICIT = true;
    if isa(varargin{2},'function_handle')
        error(nargchk(3,5,nargin));
        base = 3;
        St   = varargin{2};
    else
        base = 2;
        St   = S;
    end
    n   = varargin{base};
    
    if any(round(n) ~= n)
        error('Expecting an integer for the size');
    end

else
    n = size(S,2);
    IMPLICIT = false;
end
if nargin >= base + 1 && ~isempty( varargin{base+1} )
    tol = varargin{base+1};
end
if nargin >= base + 2 && ~isempty( varargin{base+2} )
    maxiter = varargin{base+2};
end


if length(n) == 1
    x = randn( n, 1);
else
    x = randn( n(1), n(2) );
end

cnt = 0;
nrmX = norm(x);
nrmX_old = Inf;
while abs(nrmX-nrmX_old) > tol*nrmX && cnt < maxiter
    x = x/nrmX;
    nrmX_old = nrmX;
    if ~IMPLICIT
        Sx = S*x;
    else
        Sx = S(x);
    end
    Sx = Sx / norm(Sx);
    if ~IMPLICIT
        x = S'*Sx;
    else
        x = St(Sx);
    end
    nrmX = norm(x);
    cnt = cnt+1;
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.