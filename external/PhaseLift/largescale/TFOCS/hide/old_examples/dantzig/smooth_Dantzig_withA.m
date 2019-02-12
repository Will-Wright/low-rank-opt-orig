function [ objective, gradient, dualPt ] = smooth_Dantzig_withA(A,At,DAtb,D,mu, x0, y, flag )
% [f,g] = smooth_Dantzig_withA(A,At,DAtb,D,mu,x0,y) is the smooth objective function in the
% smoothed Dantzig Selector problem, evaluated at point y.
% A,At,DATb,D,mu, and x0 are parameters
%   DATb should be D.*A'*b

% a little bit of extra overhead, but simpler for the user:
if ~isa(At,'function_handle')
    if isempty(At)
        if isa(A,'function_handle')
            error('When not providing At, A must be a matrix');
        end
        At = @(x)A'*x;
    else
        At = @(x)At*x;
    end
end
if ~isa(A,'function_handle')
    A = @(x)A*x;
end
if ~isa(D,'function_handle')
    if ~isvector(D)
        D = diag(D);
    end
    D = @(x) D.*x;
end

d = @(x) mu*norm(x-x0)^2/2;  % this is the smoothing function


AtAy = At( A( D(y) ) );

x = shrink( x0 - AtAy/mu, 1/mu ); % 
% x = shrink( x0*mu - AtAy, 1)/mu; % new, trial

if nargout == 1 && nargin == 8 && strcmpi(flag,'dual')
    % this is not the "objective", it is a dual point
    objective = x;
    return
end


objective = norm(x,1) + d(x) + (AtAy'*x-DAtb'*y);
% the objective is just the lagrangian evaluated at y and x, where x=x(y)

% my old way:
% objective = (mu/2)*( norm(x0)^2 - norm(x)^2 ) - DAtb'*y;

% MCG's new way:
% objective = (0.5/mu) * ( norm(mu_x0)^2 - norm(mu_x)^2 ) - y(end);
if nargout > 1
%     gradient = [ (1/mu) * mu_x ; -1 ];  % MCG new way
    
%     gradient = D( At( A(x) - b ) );
    gradient = D( At( A(x) )) - DAtb;
end
if nargout > 2
    dualPt = x;
end
