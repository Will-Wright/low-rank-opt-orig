rng('default')
n = 128;              % signal is n x 1 (possibly complex valued)
n_illumin = 5;       % number of structured illumination patterns

% Set signal and illumination patterns
x0      = randn(n,1) + 1i*randn(n,1);

patterntype = 'gaussian';
Pattern = MakePattern1d(n,n_illumin,patterntype);

% Generate noiseless observations
Af = @(X) FactoredIllumination(X,Pattern);

y = Af(x0(:));

%% Select model

%type = 'Poisson';
 type = 'Gaussian';

%% Select rank of matrix for truncated psd projection
RNK = 10;

%%  Set up objective function and projection routine   
lambda = 0.05; W = 1;


% Set starting point
Xinit = zeros(n,RNK);
Xinit(:,1) = randn(n,1) + 1i*randn(n,1);
Xinit = Xinit/norm(Xinit,'fro')*norm(x0(:));

% Set TFOCS options
opts.alg        = 'AT';
opts.maxIts     = 100;
opts.tol        = 1e-6;
opts.printEvery = 1;

proj = [];
epsilon = Inf;  % setting "epsilon = Inf" means W = Identity

obj = @(U) AugNegLogLikelihood(U,Af,y,type,lambda,W);

% setting additional parameters for reconstruction algorithm
pars{1} = Pattern;
pars{2} = lambda;
pars{3} = epsilon;
pars{4} = RNK;

% Call TFOCS!
tic, [X,out,opts] = tfocs(obj, proj, pars, Xinit, opts); toc
% How well we do?
fprintf('\n\nHow well do we do?\n');
rec_err = norm(x0*x0' - X(:,1)*X(:,1)')/norm(x0*x0');
fea_err = norm(Af(X(:,1))-y)/norm(y);
fprintf('relative error |x-x0| = %9.4e\n', rec_err);
fprintf('relative error |Ax-y| = %9.4e\n', fea_err);
%svd(X)'

% MPF Jan 9: Comment out for single-shot solve.
%reweightloop;




