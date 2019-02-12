function op = proj_Rplus

%PROJ_RPLUS    Projection onto the nonnegative orthant.
%    OP = PROJ_RPLUS returns an implementation of the indicator 
%    function for the nonnegative orthant.

op = @proj_Rplus_impl;

function [ v, x ] = proj_Rplus_impl( x, t )
v = 0;
switch nargin,
	case 1,
		if nargout == 2,
			error( 'This function is not differentiable.' );
		end
	case 2,
		x = max( x, 0 );
	otherwise,
		error( 'Not enough arguments.' );
end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
