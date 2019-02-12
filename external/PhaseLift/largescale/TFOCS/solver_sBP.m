function [ x,  out, opts ] = solver_sBP( A, b, mu, x0, z0, opts )

% [ x, out, opts ] = solver_sBP( A, b, mu, x0, z0, opts )
%    Solves the smoothed basis pursuit problem
%        minimize norm(x,1) + 0.5*mu*(x-x0).^2
%        s.t.     A * x == b
%    by constructing and solving the composite dual
%        maximize - g_sm(z)
%    where
%        g_sm(z) = sup_x <z,Ax-b>-norm(x,1)-(1/2)*mu*norm(x-x0)
%    A must be a linear operator or matrix, and b must be a vector. The
%    initial point x0 and the options structure opts are optional.

% Supply default values
error(nargchk(3,6,nargin));
if nargin < 4, x0 = []; end
if nargin < 5, z0 = []; end
if nargin < 6, opts = []; end
if ~isfield( opts, 'restart' ), opts.restart = 400; end


% -- legacy options from original software --
if isfield(opts,'lambda0')
    opts = rmfield(opts,'lambda0');
end
if isfield(opts,'xPlug')
    opts = rmfield(opts,'xPlug');
end
if isfield(opts,'solver')
    svr     = opts.solver;
    opts    = rmfield(opts,'solver');
    if isfield(opts,'alg') && ~isempty(opts.alg)
        disp('Warning: conflictiong options for the algorithm');
    else
        % if specified as "solver_AT", truncate:
        s = strfind( svr, '_' );
        if ~isempty(s), svr = svr(s+1:end); end
        opts.alg = svr;
    end
end




[x,out,opts] = tfocs_SCD( prox_l1, { A, -b }, proj_Rn, mu, x0, z0, opts );

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.

