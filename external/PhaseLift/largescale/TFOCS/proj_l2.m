function op = proj_l2( q )

%PROJ_LINF   The scaled 2-norm ball.
%    OP = PROJ_LINF( Q ) returns an operator implementing the 
%    indicator function for the 2-norm ball of size q,
%    { X | norm( X, 2 ) <= q }. Q is optional; if omitted,
%    Q=1 is assumed. But if Q is supplied, it must be a positive
%    real scalar.

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
op = @(varargin)proj_l2_q( q, varargin{:} );

function [ v, x ] = proj_l2_q( q, x, t )
v = 0;
switch nargin,
	case 2,
		if nargout == 2,
			error( 'This function is not differentiable.' );
		elseif norm( x, 'fro'  ) > q,
			v = Inf;
		end
	case 3,
		x = x .* ( q / norm( x, 'fro' ) );
	otherwise,
		error( 'Not enough arguments.' );
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
