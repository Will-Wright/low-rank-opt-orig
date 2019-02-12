function op = smooth_quad( P, q, r )

%SMOOTH_QUAD   Quadratic function generation.
%   FUNC = SMOOTH_QUAD( P, q, r ) returns a function handle that implements
%
%        FUNC(X) = 0.5 * TFOCS_DOT( P * x, x ) + TFOCS_DOT( q, x ) + r.
%
%   All arguments are optional; the default values are P=I, q=0, r=0. In
%   particular, calling FUNC = SMOOTH_QUAD with no arguments yields
%   
%        FUNC(X) = 0.5 * TFOCS_NORMSQ( X ) = 0.5 * TFOCS_DOT( X, X ).
%
%   If supplied, P must be a scalar, square matrix, or symmetric linear
%   operator. Furthermore, it must be positive semidefinite (convex) or
%   negative semidefinite (concave). TFOCS does not verify operator 
%   symmtetry or definiteness; that is your responsibility.

if nargin == 0,
    op = @smooth_quad_simple;
    return
end
if isa( P, 'function_handle' ),
    sz = P([],0);
    if ~isequal( sz{1}, sz{2} ),
        error( 'P must be square.' );
    end
elseif ~isnumeric( P ),
    error( 'P must be a scalar, matrix, or linear operator.' );
elseif ndims( P ) > 2 || size( P, 1 ) ~= size( P, 2 ),
    error( 'P must be a square matrix.' );
end
if nargin < 2 || isempty( q ),
    q = 0;
elseif numel(P) > 1 && numel(q) > 1 && length(q) ~= size(P,1),
    error( 'Dimension mismatch between p and q.' );
end
if nargin < 3 || isempty( r ),
    r = 0;
elseif numel(r) > 1 || ~isreal( r ),
    error( 'r must be a real scalar.' );e
end
if isnumeric( P ),
    P = 0.5 * ( P + P' );
    op = @(varargin)smooth_quad_matrix( P, q, r, varargin{:} );
else
    op = @(varargin)smooth_quad_linop( P, q, r, varargin{:} );
end

function [ v, g ] = smooth_quad_matrix( P, q, r, x, t )
if nargin == 5,
    error( 'Proximity minimization not supported by this function.' );
end
g = P * x + q;
v = 0.5 *  tfocs_dot( x, g + q ) + r;

function [ v, g ] = smooth_quad_linop( P, q, r, x, t )
if nargin == 5,
    error( 'Proximity minimization not supported by this function.' );
end
g = P( x, 1 ) + q;
v = 0.5 * tfocs_dot( x, g + q ) + r;

function [ v, x ] = smooth_quad_simple( x, t )
switch nargin,
    case 1,
    case 2,
        x = ( 1 - t ) * z;
end
v = 0.5 * tfocs_normsq( x );

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
