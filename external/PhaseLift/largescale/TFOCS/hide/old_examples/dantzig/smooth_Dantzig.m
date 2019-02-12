function [ objective, gradient, dualPt ] = smooth_Dantzig( mu, mu_x0, y, flag )
% [f,g] = smooth_Dantzig(mu,mu*x0,y) is value and gradient, resp. of the smooth objective function in the
% smoothed Dantzig Selector problem.
% [dualPt] = smooth_Dantzig(mu,mu*x0,y,flag) returns a dual point
%   (not necesssarily feasible) if "flag" is 'dual'
% [f,g,dualPt] = ...

if iscell(y)
    Ay = y{1};
    by = y{2};
else
    Ay = y(1:end-1);
    by = y(end);
end

mu_x = shrink( mu_x0 - Ay, 1 );

if nargout == 1 && nargin == 4 && strcmpi(flag,'dual')
    % this is not the "objective", it is a dual point
    objective = mu_x/mu;
    return
end

objective = (0.5/mu) * ( norm(mu_x0)^2 - norm(mu_x)^2 ) - by;

% -- An alternative -- seems to make little difference
% x=mu_x/mu; x0=mu_x0/mu;
% objective = norm(x,1)+mu/2*norm(x-x0)^2+y(1:end-1)'*x-y(end) ;
if nargout > 1,
    if iscell(y)
        gradient = (1/mu)*mu_x;
    else
        gradient = [ (1/mu) * mu_x ; -1 ];
    end
end

if nargout > 2
    dualPt = mu_x/mu;
end
