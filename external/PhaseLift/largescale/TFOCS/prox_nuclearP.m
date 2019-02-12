function OP = prox_nuclearP( q )

%PROX_NUCLEARP    Nuclear norm.
%    OP = PROX_NUCLEARP( q ) implements the nonsmooth function
%        OP(X) = q * sum(svd(X)).
%    Q is optional; if omitted, Q=1 is assumed. But if Q is supplied, 
%    it must be a positive real scalar.
%
% This implementation uses the packSVD object so that the low rank
% structure of the iterates may be preserved.

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
%op = @(x,t)prox_nuclearP_impl( q, x, t );
op = @(varargin)prox_nuclearP_impl( q, varargin{:} );

function [ v, X ] = prox_nuclearP_impl( q, X, t )
if nargin > 2 && t > 0,
    X = svd_shrink( X, q * t );
end
v = q * sum(svd(X));

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
