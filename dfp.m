function [y,stats] = dfp(A, b, e, s, x, y0, varargin)
%DFP  Dual from primal.
import util.*

p = inputParser;
p.addParameter('gTol', 1e-6);
p.addParameter('fTol', -inf);
p.addParameter('verbose', false);
p.addParameter('maxit', 1e3);
p.parse(varargin{:});

l    = s/dot(x(:),x(:));
v    = x/norm(x(:));
xlhs =   v;
xrhs = l*v;

if (nargin <= 5) || isempty(y0)
   y0 = project(zeros(size(b)),b,e,s);
end

P  = @(Y)project(Y,b,e,s);
% FG = @(Y)A.dfpobjective(xlhs,xrhs,Y);
FG = A.dfpobjective(xlhs,xrhs);

pgtol   = p.Results.gTol;
fTol    = p.Results.fTol;
verbose = p.Results.verbose;
maxit   = p.Results.maxit;

% Tol for small residual. SPG objective is 1/4||r||^2, so adjust.
if fTol > 0
   ftol = fTol*(1+normv(xrhs));
   ftol = ftol^2/2;
else
   ftol = fTol;
end

[y,stats] = spg(FG,P,y0,struct('verbose',verbose,'pgtol',pgtol,...
                'ftol',ftol,'maxit',maxit));
end
