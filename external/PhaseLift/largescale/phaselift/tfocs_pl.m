function [x, r, info] = tfocs_pl(A, b, varargin)
%TFOCS_PL  Solve the PhaseLift problem using TFOCS.
%
% [x, r, info] = tfocs_pl(A, b)  where A is a PL operator.
%
% Depends on code in external/PhaseLift/largescale/phaselift/.

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('fid', 1);
ip.addParameter('verbosity', 1);
ip.addParameter('iterations', 1000);
ip.parse(varargin{:});

n = A.n;
L = A.m3;
AdjGradient();         % reset the adjoint-call counter
FactoredIllumination();% reset the adjoint-call counter
M = squeeze(A.masks);  % TFOCs will want a matrix for the 1D case
y = reshape(b, n, L);  % TFOCS needs different size RHS
RNK = 10;              % rank of matrix for truncated psd projection
lambda = 0.05; W = 1;  % objective function and projection routine

% Set starting point
Xinit = zeros(n,RNK);
Xinit(:,1) = randn(n,1) + 1i*randn(n,1);
Xinit = Xinit/norm(Xinit,'fro');

% Set TFOCS options
opts.fid        = ip.Results.fid;
opts.alg        = 'AT';
opts.maxIts     = ip.Results.iterations;
opts.tol        = 0;%1e-6;
opts.printEvery = ip.Results.verbosity;
opts.saveHist   = false;
opts.countOps = true;

proj = [];
epsilon = Inf;  % setting "epsilon = Inf" means W = Identity

Af = @(X) FactoredIllumination(X, M);
obj = @(U) AugNegLogLikelihood(U, Af, y, 'gaussian', lambda, W);

% setting additional parameters for reconstruction algorithm
pars{1} = M;
pars{2} = lambda;
pars{3} = epsilon;
pars{4} = RNK;

% Call TFOCS!
tStart = tic;
[X, info] = tfocs(obj, proj, pars, Xinit, opts);

% Compute output quantities.
x = X(:,1); %[U,S,V] = svd(X); x = S(1,1)*U(:,1);
r = b - A.forward(x);
info.time = toc(tStart);
info.nAdjoint = AdjGradient;
info.nMeasure = RNK*FactoredIllumination;
info.nfft = 2*L*info.nAdjoint + ...  % 2L ffts per adjoint
            RNK*L*info.nMeasure;     %  L ffts per RNK per forward 

end
