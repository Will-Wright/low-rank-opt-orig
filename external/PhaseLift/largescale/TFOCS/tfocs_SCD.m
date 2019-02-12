function [ x, odata, opts ] = tfocs_SCD( objectiveF, affineF, dualproxF, mu, x0, z0, opts )

% [ x, out, opts ] = tfocs_SCD( objectiveF, affineF, dualproxF, x0, z0, opts )
%   Solves a conic problem using the smoothed conic dual approach. The goal
%   is to solve a problem in the following conic form:
%       minimize objectiveF(x)+0.5*mu(x-x0).^2
%       s.t.     affineF(x) \in \cK
%   The user is responsible for constructing the dual proximity function
%   so that the dual can be described by the saddle point problem
%        maximize_z inf_x [ objectiveF(x)+0.5*mu(x-x0).^2-<affineF(x),z> ] -dualproxF(z)
%   If mu = Inf, the method ignores the objective function and solves
%       minimize 0.5*(x-x0).^2
%       s.t.     affineF(x) \in \cK

if nargin == 0
    if nargout > 0, x = tfocs();
    else tfocs();
    end
    return;
end

error(nargchk(4,7,nargin));
if nargin < 5, x0 = []; end
if nargin < 6, z0 = []; end
if nargin < 7, opts = []; end

% Handle special cases of zero objective of infinite mu
if isinf( mu ),
    mu = 1;
    objectiveF = prox_0;
elseif isempty( objectiveF ),
    objectiveF = prox_0;
elseif iscell( objectiveF )  % allow the case of {[],[],...,[]}
    for k = 1:length(objectiveF)
        if isempty( objectiveF{k} ) || isscalar(objectiveF{k})
            objectiveF{k} = prox_0;
        end
    end
end

% The affine quantities will be used in adjoint orientation
if isfield( opts, 'adjoint'  ),
    opts.adjoint = ~opts.adjoint;
else
    opts.adjoint = true;
end
opts.saddle = true;
opts.maxmin = -1;
if isempty( x0 ), x0 = 0; end
if iscell(objectiveF)
    for k = 1 : length(objectiveF)
        smoothF{k} = @(varargin)smooth_dual( objectiveF{k}, 1/mu, x0, varargin{:} );
    end
else
    smoothF = @(varargin)smooth_dual( objectiveF, 1/mu, x0, varargin{:} );
end

if isempty(dualproxF)
    dualproxF = proj_Rn;
elseif iscell(dualproxF)
    for k = 1:length(dualproxF)
        dp = dualproxF{k};
        if isempty( dp ) || isnumeric( dp ) && numel( dp ) == 1 && dp == 0,
            dualproxF{k} = proj_Rn;
        end
    end
end
try  % this is annoying for debugging
    [ z, odata, opts ] = tfocs( smoothF, affineF, dualproxF, z0, opts );
catch err
    if strfind(err.message,'x0')
        fprintf(2,'Error involves z0 (which is referred to as x0 below)\n');
    end
    rethrow(err);
end
opts.adjoint = ~opts.adjoint;
opts = rmfield( opts, { 'saddle', 'maxmin' } );
x = odata.dual;
odata.dual = z;

function [ prox, x ] = smooth_dual( objectiveF, mu_i, x0, ATz )
% Adding 0 to ATz will destroy the sparsity
if (isscalar(x0) && x0 == 0) || numel(x0) == 0 || nnz(x0) == 0
    [ v, x ] = objectiveF( mu_i * ATz, mu_i );
else
    [ v, x ] = objectiveF( x0 + mu_i * ATz, mu_i );
end
prox = tfocs_dot( ATz, x ) - v - (0.5/mu_i) * tfocs_normsq( x - x0 );
prox = -prox;
x = -x;

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.



