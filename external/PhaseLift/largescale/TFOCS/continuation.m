function [ x, odata, optsOut ] = continuation( fcn, mu, x0, z0, opts, contOpts )
% [...] = CONTINUATION( FCN, MU, X0, Z0, OPTS, CONT_OPTS )
%   is a wrapper to perform continuation on FCN, where
%   FCN is a function that calls tfocs_SCD or a tfocs solver.
%   FCN must accept input arguments MU, X0, Z0, and OPTS,
%       e.g. FCN = @(mu,x0,z0,opts) solver_sBP( A, b, mu, x0, z0, opts )
%
%   CONT_OPTS are options that affect how continuation is performed.
%   To see the default options, call this script with no inputs.
%
%   Options for CONT_OPTS:
%       maxIts      - max # of continuation iterations
%       accel       - use accelerated continuation or not
%       betaTol     - every continuation iteration, the tolerance
%                       is decreased by 'betaTol'
%       innerTol    - every continuation iteration, except the last one,
%                       is solved to this tolerance.
%                       This option overrides 'betaTol'
%       tol         - the outer loop will stop when either
%                       'maxIts' is reached, or when the change
%                       from one step to the next is less than 'tol'.
%       innerMaxIts - maximum number of inner iterations during
%                       every continuation iteration except the last one.

if nargin < 6, contOpts = []; end
kMax    = setOpts('maxIts',3);
ACCEL   = setOpts('accel',true);
betaTol = setOpts('betaTol',2);
innerTol= setOpts('innerTol',[] );
stopTol = setOpts('tol',1e-3);
innerMaxIts = setOpts('innerMaxIts',[] );

if nargin == 0
    if nargout == 0
        disp('Default options for CONTINUATION are:');
    end
    x = contOpts;
    return;
end

error(nargchk(2,6,nargin));
if nargin < 3, x0   = []; end
if nargin < 4, z0   = []; end
if nargin < 5, opts = []; end

% what is the user specified stopping tolerance?
if ~isfield(opts,'tol')
    % get the default value:
    defaultOpts = tfocs_SCD();
    opts.tol = defaultOpts.tol;
end
if isfield(opts,'fid'), fid = opts.fid; else fid = 1; end
finalTol    = opts.tol;
tol         = finalTol*betaTol^(kMax+1);

xOld = x0;
odataCumulative      = [];
odataCumulative.niter= 0;
for k = 1:kMax
    if kMax > 1
        fprintf(fid,'---- Continuation step %d of %d ----\n', k, kMax );
    end
    
    optsTemp    = opts;
    if ~isempty( innerMaxIts ) && k < kMax
        optsTemp.maxIts     = innerMaxIts;
    end
    tol             = tol/betaTol;
    if ~isempty( innerTol ) && k < kMax
        optsTemp.tol    = innerTol;
    else
        optsTemp.tol    = tol;
    end
    
    % call the solver
    [x, odata, optsOut ] = fcn( mu, x0, z0, optsTemp );
    
    % update output data
    fields = { 'f', 'normGrad', 'stepsize','theta','counts','err' };
    for ff = fields
        f = ff{1};
        if isfield(odata,f) 
            if isfield(odataCumulative,f)
                odata.(f) = [odataCumulative.(f); odata.(f)];
            end
            odataCumulative.(f) = odata.(f);
        end
        
    end
    if isfield(odata,'niter')
        % include the old iterations in the count
        odataCumulative.niter = odata.niter + odataCumulative.niter;
        odata.niter     = odataCumulative.niter;
        
        % record at which iterations continuation happened
        if isfield(odataCumulative,'contLocations')
            odataCumulative.contLocations = [odataCumulative.contLocations; odata.niter];
            odata.contLocations = odataCumulative.contLocations;
        else
            odata.contLocations = odata.niter;
            odataCumulative.contLocations = odata.niter;
        end
    end
    
    if k == kMax, break; end
    
    % Update the prox center
    if isempty(xOld), xOld = zeros(size(x)); end
    if ACCEL
        x0 = x + (k-1)/(k+2)*( x - xOld );
    else
        x0 = x;
    end
    
    if isa( odata.dual, 'tfocs_tuple')
        z0 = cell( odata.dual );
    else
        z0 = odata.dual;
    end

    if norm(x - xOld)/norm(xOld) <= stopTol
        fprintf(fid,'Continuation ending, due to convergence\n');
        break;
    end
    xOld = x;
    if k == kMax
        fprintf(fid,'Continuation ending, due to reaching maximum number of outer iterations\n');
    end
end


function out = setOpts(fieldName,default)
    if ~isfield(contOpts,fieldName)
        contOpts.(fieldName) = default;
    end
    out = contOpts.(fieldName);
end


end

% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.



