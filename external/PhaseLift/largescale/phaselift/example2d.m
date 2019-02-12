n = 32;              % signal is n x n (possibly complex valued)
n_illumin = 8;       % number of structured illumination patterns

% Set signal and illumination patterns
x0      = randn(n,n) + 1i*randn(n,n);


Pattern = randn(n,n,n_illumin) + 1i*randn(n,n,n_illumin);
%Pattern = round(rand(n,n,n_illumin));

% Generate noiseless observations
Af = @(X) FactoredIllumination(X,Pattern);

y = Af(x0(:));

SNR = 0.05;
ynoise = abs(randn(size(y)));
ynoise = ynoise/normfro(ynoise)*normfro(y)*SNR;

y = y + ynoise;

%% Select model

% type = 'Poisson';
 type = 'Gaussian';  % although noise is not truely Gaussian 

%% Select rank of matrix for truncated psd projection
RNK = 10;

%%  Set up objective function and projection routine   
lambda = 0.05; W = 1;


% Set starting point
%Xinit = randn(n^2,RNK) + i*randn(n^2,RNK);
%Xinit = Xinit/norm(Xinit,'fro')*norm(x0(:));
Xinit1 = randn(n^2,1); Xinit2 = Xinit1(1:RNK); Xinit = Xinit1*Xinit2';
Xinit = Xinit/norm(Xinit(:))*norm(x0(:));



% Set TFOCS options
opts.alg        = 'AT';
opts.maxIts     = 200;
opts.tol        = 1e-6;
opts.printEvery = 50;
opts.maxmin     = 1; % Maximize

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
fprintf('Approximate relative error = %9.4e\n', phaserr(x0(:),X(:,1)));
svd(X)'
norm(Af(X(:,1))-y)/norm(y)


%reweightloop;




