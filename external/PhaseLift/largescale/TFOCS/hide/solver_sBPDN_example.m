function [ x, z, odata, opts ] = solver_sBPDN_example( A, b, epsilon, mu, x0, z0, opts )

% [ x, z, out, opts ] = solver_sBPDN_example( A, b, epsilon, mu, x0, opts )
%    Solves the smoothed basis pursuit denoising problem
%        minimize norm(x,1) + 0.5*mu*(x-x0).^2
%        s.t.     norm(A*x-b,2) <= epsilon
%    by constructing and solving the composite dual
%        maximize - g_sm(z) - epsilon*norm(z,2)
%    where
%        gsm(z) = sup_x <z,Ax-b>-norm(x,1)-(1/2)*mu*norm(x-x0)
%    A must be a linear operator or matrix, and b must be a vector. The
%    initial point x0 and the options structure opts are optional.
%
%   This file is similar to solver_sBPDN.m, except that instead of calling
%   tfocs_SCD, it calls tfocs directly.  This is to make it more transparent,
%   and serve as a template for more complicated problems.

% Supply default values
error(nargchk(4,7,nargin));
if nargin < 4, x0 = []; end
if nargin < 5, z0 = []; end
if nargin < 7, opts = []; end
if ~isfield( opts, 'restart' ), opts.restart = 400; end

% "A" should be a matrix or function handle that conforms to SPARCO conventions
% If the user has given instead {Af,At,M,N}, then we will put it in SPARCO-compatible
% form for them
if iscell(A)
    [M,N,Af,At] = deal( A );
    if ~isscalar(M) || ~isscalar(N) || ~isa(Af,'function_handle') || ~isa(At,'function_handle')
        error('If "A" is a cell array, it should be in format {M,N,Af,At}');
    end
    %A = @(x,mod)linop_SPARCO(M,N,Af,At,x,mode);
    A = linop_handles([M,N],Af,At);
end

%[x,out,opts] = tfocs_SCD( @prox_l1, { A, -b }, prox_scale( @prox_l2, epsilon ), mu, x0, z0, opts );

objectiveF  = prox_l1;
affineF     = {A,-b};
if epsilon
    dualproxF   = prox_l2( epsilon );
else
    dualproxF   = proj_Rn;
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
smoothF = @(varargin)smooth_dual( objectiveF, mu, x0, varargin{:} );

% Call the solver
[ z, odata, opts ] = tfocs( smoothF, affineF, dualproxF, z0, opts );

opts.adjoint = ~opts.adjoint;
opts = rmfield( opts, { 'saddle', 'maxmin' } );
x = odata.dual;
odata.dual = z;

% if objectiveF = @prox_l1, this is what we do:
function [ prox, x ] = smooth_dual( objectiveF, mu, x0, ATz )
[ fx, x ] = objectiveF( x0 + (1/mu)*ATz, mu );
prox = -tfocs_dot( ATz, x ) + fx + 0.5 * mu * tfocs_normsq( x - x0 );
x = -x;



% Copyright 2010 Michael C. Grant and Stephen R. Becker.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

