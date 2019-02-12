function op = prox_l2( q )

%PROX_L2    L2 norm.
%    OP = PROX_L2( q ) implements the nonsmooth function
%        OP(X) = q * norm(X,'fro').
%    Q is optional; if omitted, Q=1 is assumed. But if Q is supplied,
%    then it must be a positive real scalar.

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
op = @(varargin)prox_l2_q( q, varargin{:} );

function [ v, x ] = prox_l2_q( q, x, t )
if nargin < 2,
	error( 'Not enough arguments.' );
end
v = sqrt( tfocs_normsq( x ) );
if nargin == 3,
	s = 1 - 1 ./ max( v / ( t * q ), 1 );
	x = x * s;
	v = v * s;
elseif nargout == 2,
	error( 'This function is not differentiable.' );
end
v = q * v;

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
