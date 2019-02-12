function [x, stats] = pfd(A, b, e, y, v, g, varargin)
%PFD  Primal from dual.
import util.*

p = inputParser;
p.addParameter('gTol', 1e-8);
p.addParameter('fTol', 1e-6);
p.addParameter('verbose', false);
p.addParameter('maxit', 1e3);
p.addParameter('solver','spg');
p.parse(varargin{:});

gTol = p.Results.gTol;
fTol = p.Results.fTol;
verb = p.Results.verbose;
maxit = p.Results.maxit;
solver = p.Results.solver;

% Parse options.
if ~isempty(e)  &&  e > 0
   brhs = b-e*y/normv(y);
else
   brhs = b;
end
if (nargin <= 5) || isempty(v) || isempty(g)
   [~,g,V,~] = A.gdobjective(y);
   v  = V(:,1);
end
% % nx = realsqrt(A.n*sum(brhs(:)))/norm(A.masks(:)); % WFlow paper's version
nx = sqrt(max(0, rdot(g,brhs))) / normv(g);
x0 = v*nx;
FG = A.pfdobjective(brhs);

% Tol for small residual. SPG objective is 1/4||r||^2, so adjust.
if fTol > 0
   fTolr = fTol*(1+normv(brhs));
   fTolr = (fTolr/2)^2; % == ||r|| <= fTolr*(1 + ||b||)
else
   fTolr = fTol;
end

[x, stats] = subsolver(FG, x0, solver, verb, gTol, fTolr, maxit);

end % function

function [x, stats] = subsolver(FG, x0, solver, verb, gTol, fTol, maxit)

if strcmp(solver, 'spg')
   [x, stats] = spg(FG,[],x0,...
      struct('verbose',verb,'pgtol',gTol,'ftol',fTol,'maxit',maxit));
   
else
   
   opts.optTol = gTol;
   opts.progTol = fTol;
   opts.MaxIter = maxit;
   opts.MaxFunEvals = inf;
   if verb
      opts.Display = 'iter';
   else
      opts.Display = 'off';
   end
   
   [ri, ~, ~, output] = minFunc( @(ri)mfPFDobjective(FG,ri), c2r(x0), opts);
   x = r2c(ri);
   [f,g] = FG(x);
   x = reshape(x,size(g));
   stats.fOut = f;
   stats.pgNorms = util.normv(g);
   stats.iterations = output.iterations;
   
end

end

function ri = c2r(z)
ri = [real(z(:));imag(z(:))];
end

function z = r2c(ri)
n = numel(ri);
z = complex(ri(1:n/2),ri(n/2+1:end));
end

function [f,g] = mfPFDobjective(FG,ri)
z = r2c(ri);
if nargout == 1
   f = FG(z);
else
   [f,g] = FG(z);
   g = c2r(g);
end
end