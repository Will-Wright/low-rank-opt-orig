function [ x, z, out, opts ] = solver_sL1Analysis( W, A, b, epsilon, mu, x0, z0, opts )

%[ x, opts ] = solver_sDantzig( A, b, delta, mu, x0, opts )
%    Solves the smoothed Dantzig
%        minimize norm(W*x,1) + (1/2)*mu*norm(x-x0).^2
%        s.t.     norm(A*x-b,2) <= epsilon
%    using the Auslender/Teboulle variant applied to the dual problem. The
%    corresponding composite dual problem is
%        maximize - g_sm(z1,n2) - epsilon*norm(z2,1)
%        s.t      norm(z1,Inf) <= 1
%    where
%        gsm(z) = sup_x <z1,Wx>+<z2,A*x-b>-(1/2)*mu*norm(x-x0)
%    A must be a linear operator, b must be a vector, and delta and mu
%    must be positive scalars. Initial points x0 and z0 are optional.
%    The standard calling sequence assumes that D=I. To supply a scaling,
%    pass the cell array { D, A } instead of A. Note that D must be a
%    vector, not a matrix or operator.

% Supply default values
error(nargchk(4,7,nargin));
if nargin < 5, x0 = []; end
if nargin < 6, z0 = []; end
if nargin < 7, opts = []; end
if ~isfield( opts, 'alg' ), opts.alg = 'AT'; end

% Build the composite linear operator and projection
WA = { W, 0 ; A, -b };
PP = { proj_linf( 1 ) ; prox_l2( epsilon ) };
[x,z,out,opts] = tfocs_SCD( [], WA, PP, mu, x0, z0, opts );

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

