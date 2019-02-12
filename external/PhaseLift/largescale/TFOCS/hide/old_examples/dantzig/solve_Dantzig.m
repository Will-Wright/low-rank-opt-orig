function [x,lambda,out,actualOpts]=solve_Dantzig( A, b, D, mu0, EPS, opts )
% x = solve_Dantzig( A, b, D, mu, epsilon )
%   solves the Basis Pursuit or Basis Pursuit De-Noising problem:
%       min_x ||x||_1
%   subject to
%       || DA'(Ax-b) ||_\infty <= EPS
%
%   A is either a matrix or a 2-entry cell array, where the 1st entry
%       is a function handle to compute A*x, and the 2nd entry is a
%       function handle to compute A'*x.
%
%   D is a scaling vector (or diagonal matrix)
%
%   epsilon is 0 by default
%
%   mu is the smoothing parameter used by the conic dual solvers
%
%   x = solve_BPDN( A, b, mu, epsilon, opts )
%   uses information in the structure "opts" to pass to the solver
%
%   [x,lambda,out] = solve_BPDN(...) returns additional information
%
%   TFOCS ver
%
%   See also solver_AT solve_BPDN
if nargin == 0
    % output the possible options
    disp('Possible fields for the "opts" structure are:');
    disp('(note: this wrapper uses different default values than those listed below)');
    solver_AT()
    return;
end
nargchk( 5,6,nargin );

mu = mu0;
if nargin < 6
    opts = [];
end

% A is a M x N matrix
if iscell(A),
    Af = A{1};  % forward multiply
    At = A{2};  % transpose multiply
    % need to know what N is:
    if isfield(opts,'N') && isnumeric(opts.N) && opts.N > 0
        N = opts.N;
    elseif ~isfield(opts,'x0')
        N = length( At(b) );
    end
    M = length(b);
else
    Af = @(x) A*x;
    At = @(x) A'*x;
    [M,N] = size(A);
end
if isfield(opts,'x0')
    xPlug = opts.x0;
    N = length(xPlug);
elseif isfield(opts,'xPlug')
    xPlug = opts.xPlug;
    N = length(xPlug);
else
    xPlug = zeros(N,1);
end

if isempty(D)
    D = ones(N,1);
elseif ~isvector(D)
    D = diag(D);
end
delta0 = EPS;

% obj    = @(x) norm(x,1) + (mu/2)*norm(x-xPlug)^2;
linear = @(varargin) linear_Dantzig( Af, At, b, D, varargin{:} );
smooth = @(varargin) smooth_Dantzig( mu, mu * xPlug, varargin{:} ); % mu*xPlug??? Yes.
proj   = @(varargin) nonsmooth_Dantzig(delta0, varargin{:} );

updateOpts( 'maxcounts', [Inf,Inf,Inf,Inf] );
updateOpts( 'maxits', 3000 );
updateOpts( 'maxmin', -1 );
updateOpts( 'saddle', 1 );
updateOpts( 'printEvery', 500 );
updateOpts( 'tol', 1e-4 );
if ~iscell(A)
%     updateOpts( 'Lexact', norm(A*A')/mu ); % for LASSO
    updateOpts( 'Lexact', norm(diag(D)*A'*A)^2/mu ); % for Dantzig Selector
end



% solver = @solver_AT;
% solver = @solver_N83;
% solver = @solver_N07;

updateOpts( 'lambda0', zeros(N,1) );

updateOpts('solver','solver_AT');  % default: use first-order solver

if strcmpi( opts.solver, 'cvx') && ~iscell(A)
    [lambda,x] = cvx_dantzig_solver(A,D,b,xPlug,delta0,N,mu);
    out = [];
else    
    solver = str2func( opts.solver );
    
    updateOpts( 'restart', Inf );
    if isempty(opts.restart), opts.restart = Inf; end
    
    maxits = opts.maxits;
    
    restart = min(opts.restart,maxits);
    
    opts.maxits = restart;
    
    lambda = opts.lambda0;
%     out = struct([]);
    out = struct('niter',0);  % weird bug requires struct to be non-empty
    for k = 1:round(maxits/restart)
        % only lambda0 is updated -- no change to x
        [lambda,out_i,actualOpts] = solver( smooth, linear, proj, lambda, opts );
        
        % copy "out" to a combined structure
        oldIter = out.niter;
        niter = out_i.niter;
        for f = fieldnames(out_i)'
            ff = f{1};
            if size( out_i.(ff), 1 ) == niter
                % for these type of fields, append
                if isfield(out,ff)
                    out.(ff) = [out.(ff); out_i.(ff)];
                else
                    out.(ff) = out_i.(ff);
                end
            else
                % for these fields, update and overwrite
                out.(ff) = out_i.(ff);
            end 
        end
        % and this special field needs to add:
        out.niter = oldIter + niter;
        
    end
    actualOpts.restart = restart;
    actualOpts.maxits = maxits;
    
    
    x = out.dual;
end



% Internal function (this has same scope as parent function)
function updateOpts(field,value)
% Warning: this is a simple version of this function,
%   it doesn't check for capitalization mismatches!
  if isfield(opts,field) && ~isempty(opts.(field))
      % do nothing: the user has already supplied a value
  else
      opts.(field) = value;
  end
end


end

% External functions:
function [lambda,x] = cvx_dantzig_solver(A,D,b,xPlug,delta0,N,mu)
if delta0    
    cvx_begin
          cvx_precision best
%            cvx_quiet true
           variable x(N,1)
           dual variable lambda
           minimize norm(x,1) + mu/2*sum_square(x-xPlug)
           subject to
                abs(diag(D)*A'*(b-A*x)) <= delta0 : lambda
    cvx_end
else 
    cvx_begin
          cvx_precision best
%            cvx_quiet true
           variable x(N,1)
           dual variable lambda
           minimize norm( x, 1 ) + mu/2*sum_square(x-xPlug)
           subject to
                b == A*x : lambda
    cvx_end
end

end
